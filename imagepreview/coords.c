/*** File libwcs/hget.c
 *** August 22, 2007
 *** By Doug Mink, dmink@cfa.harvard.edu
 *** Harvard-Smithsonian Center for Astrophysics
 *** Copyright (C) 1994-2007
 *** Smithsonian Astrophysical Observatory, Cambridge, MA, USA
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
 Correspondence concerning WCSTools should be addressed as follows:
 Internet email: dmink@cfa.harvard.edu
 Postal address: Doug Mink
 Smithsonian Astrophysical Observatory
 60 Garden St.
 Cambridge, MA 02138 USA
 
 * Module:	hget.c (Get FITS Header parameter values)
 * Purpose:	Extract values for variables from FITS header string
 * Subroutine:	hgeti2 (hstring,keyword,ival) returns short integer
 * Subroutine:	hgeti4c (hstring,keyword,wchar,ival) returns long integer
 * Subroutine:	hgeti4 (hstring,keyword,ival) returns long integer
 * Subroutine:	hgetr4 (hstring,keyword,rval) returns real
 * Subroutine:	hgetra (hstring,keyword,ra) returns double RA in degrees
 * Subroutine:	hgetdec (hstring,keyword,dec) returns double Dec in degrees
 * Subroutine:	hgetr8c (hstring,keyword,wchar,dval) returns double
 * Subroutine:	hgetr8 (hstring,keyword,dval) returns double
 * Subroutine:	hgetl  (hstring,keyword,lval) returns logical int (0=F, 1=T)
 * Subroutine:	hgetsc (hstring,keyword,wchar,lstr,str) returns character string
 * Subroutine:	hgets  (hstring,keyword, lstr, str) returns character string
 * Subroutine:	hgetm  (hstring,keyword, lstr, str) returns multi-keyword string
 * Subroutine:	hgetdate (hstring,keyword,date) returns date as fractional year
 * Subroutine:  hgetndec (hstring, keyword, ndec) returns number of dec. places
 * Subroutine:	hgetc  (hstring,keyword) returns character string
 * Subroutine:	blsearch (hstring,keyword) returns pointer to blank lines
 before keyword
 * Subroutine:	ksearch (hstring,keyword) returns pointer to header string entry
 * Subroutine:	str2ra (in) converts string to right ascension in degrees
 * Subroutine:	str2dec (in) converts string to declination in degrees
 * Subroutine:	strsrch (s1, s2) finds string s2 in null-terminated string s1
 * Subroutine:	strnsrch (s1, s2, ls1) finds string s2 in ls1-byte string s1
 * Subroutine:	hlength (header,lhead) sets length of FITS header for searching
 * Subroutine:  isnum (string) returns 1 if integer, 2 if fp number, else 0
 * Subroutine:  notnum (string) returns 0 if number, else 1
 * Subroutine:  numdec (string) returns number of decimal places in numeric string
 * Subroutine:	strfix (string,blankfill,zerodrop) removes extraneous characters
 */

#include <string.h>		/* NULL, strlen, strstr, strcpy */
#include <stdio.h>
#include <stdlib.h>
#ifndef VMS
#include <limits.h>
#else
#define INT_MAX  2147483647 /* Biggest number that can fit in long */
#define SHRT_MAX 32767
#endif
#define VLENGTH 81

static char val[VLENGTH+1];
static int multiline = 0;

double str2ra (const char *in);
double str2dec (const char *in);
char *strsrch (const char *s1, const char *s2);
char *strnsrch (const char *s1, const char *s2, const int ls1);
char *strncsrch (const char *s1, const char *s2, const int ls1);


/* Return the right ascension in degrees from sexagesimal hours or decimal degrees */

double str2ra (const char *in)
{
    double ra;	/* Right ascension in degrees (returned) */
    
    ra = str2dec (in);
    if (strsrch (in,":"))
        ra = ra * 15.0;
    
    return (ra);
}


/* Return the declination in degrees from sexagesimal or decimal degrees */

double str2dec (const char *in)
{
    double dec;		/* Declination in degrees (returned) */
    double deg, min, sec, sign;
    char *value, *c1, *c2;
    int lval;
    char *dchar;
    
    dec = 0.0;
    
    /* Return 0.0 if string is null */
    if (in == NULL)
        return (dec);
    
    /* Translate value from ASCII colon-delimited string to binary */
    if (in[0]) {
        value = (char *) in;
        
        /* Remove leading spaces */
        while (*value == ' ')
            value++;
        
        /* Save sign */
        if (*value == '-') {
            sign = -1.0;
            value++;
	    }
        else if (*value == '+') {
            sign = 1.0;
            value++;
	    }
        else
            sign = 1.0;
        
        /* Remove trailing spaces */
        lval = strlen (value);
        while (value[lval-1] == ' ')
            lval--;
        
        if ((c1 = strsrch (value,":")) == NULL)
            c1 = strnsrch (value," ",lval);
		
        if (c1 != NULL) {
            *c1 = 0;
            
            deg = (double) atoi (value);
            *c1 = ':';
            value = c1 + 1;
            if ((c2 = strsrch (value,":")) == NULL)
                c2 = strsrch (value," ");
            if (c2 != NULL) {
                *c2 = 0;
                min = (double) atoi (value);
                *c2 = ':';
                value = c2 + 1;
                sec = atof (value);
                
            }
            else {
                sec = 0.0;
                if ((c1 = strsrch (value,".")) != NULL)
                    min = atof (value);
                if (strlen (value) > 0)
                    min = (double) atoi (value);
            }
            dec = sign * (deg + (min / 60.0) + (sec / 3600.0));
	    }
        else if (isnum (value) == 2) {
            if ((dchar = strchr (value, 'D')))
                *dchar = 'e';
            if ((dchar = strchr (value, 'd')))
                *dchar = 'e';
            if ((dchar = strchr (value, 'E')))
                *dchar = 'e';
            dec = sign * atof (value);
	    }
        else
            dec = sign * (double) atoi (value);
	}
    return (dec);
}


/* Find string s2 within null-terminated string s1 */

char *
strsrch (s1, s2)

const char *s1;	/* String to search */
const char *s2;	/* String to look for */

{
    int ls1;
    ls1 = strlen (s1);
    return (strnsrch (s1, s2, ls1));
}


/* Find string s2 within string s1 */

char *
strnsrch (s1, s2, ls1)

const char *s1;	/* String to search */
const char *s2;	/* String to look for */
const int ls1;	/* Length of string being searched */

{
    char *s,*s1e;
    char cfirst,clast;
    int i,ls2;
    
    /* Return null string if either pointer is NULL */
    if (s1 == NULL || s2 == NULL)
        return (NULL);
    
    /* A zero-length pattern is found in any string */
    ls2 = strlen (s2);
    if (ls2 ==0)
        return ((char *) s1);
    
    /* Only a zero-length string can be found in a zero-length string */
    if (ls1 ==0)
        return (NULL);
    
    cfirst = (char) s2[0];
    clast = (char) s2[ls2-1];
    s1e = (char *) s1 + (int) ls1 - ls2 + 1;
    s = (char *) s1;
    while (s < s1e) {
        
        /* Search for first character in pattern string */
        if (*s == cfirst) {
            
            /* If single character search, return */
            if (ls2 == 1)
                return (s);
            
            /* Search for last character in pattern string if first found */
            if (s[ls2-1] == clast) {
                
                /* If two-character search, return */
                if (ls2 == 2)
                    return (s);
                
                /* If 3 or more characters, check for rest of search string */
                i = 1;
                while (i < ls2 && s[i] == s2[i])
                    i++;
                
                /* If entire string matches, return */
                if (i >= ls2)
                    return (s);
            }
	    }
        s++;
	}
    return (NULL);
}


/* Find string s2 within null-terminated string s1 (case-free search) */

char *
strcsrch (s1, s2)

const char *s1;	/* String to search */
const char *s2;	/* String to look for */

{
    int ls1;
    ls1 = strlen ((char *) s1);
    return (strncsrch (s1, s2, ls1));
}


/* Find string s2 within string s1 (case-free search) */

char *
strncsrch (s1, s2, ls1)

const char *s1;	/* String to search */
const char *s2;	/* String to look for */
const int ls1;	/* Length of string being searched */

{
    char *s,*s1e, sl, *os2;
    char cfirst,ocfirst;
    char clast = ' ';
    char oclast = ' ';
    int i,ls2;
    
    /* Return null string if either pointer is NULL */
    if (s1 == NULL || s2 == NULL)
        return (NULL);
    
    /* A zero-length pattern is found in any string */
    ls2 = strlen (s2);
    if (ls2 ==0)
        return ((char *) s1);
    
    /* Only a zero-length string can be found in a zero-length string */
    os2 = NULL;
    if (ls1 ==0)
        return (NULL);
    
    /* For one or two characters, set opposite case first and last letters */
    if (ls2 < 3) {
        cfirst = (char) s2[0];
        if (cfirst > 96 && cfirst < 123)
            ocfirst = cfirst - 32;
        else if (cfirst > 64 && cfirst < 91)
            ocfirst = cfirst + 32;
        else
            ocfirst = cfirst;
        if (ls2 > 1) {
            clast = s2[1];
            if (clast > 96 && clast < 123)
                oclast = clast - 32;
            else if (clast > 64 && clast < 91)
                oclast = clast + 32;
            else
                oclast = clast;
	    }
	}
    
    /* Else duplicate string with opposite case letters for comparison */
    else {
        os2 = (char *) calloc (ls2, 1);
        for (i = 0; i < ls2; i++) {
            if (s2[i] > 96 && s2[i] < 123)
                os2[i] = s2[i] - 32;
            else if (s2[i] > 64 && s2[i] < 91)
                os2[i] = s2[i] + 32;
            else
                os2[i] = s2[i];
	    }
        cfirst = s2[0];
        ocfirst = os2[0];
        clast = s2[ls2-1];
        oclast = os2[ls2-1];
	}
    
    /* Loop through input string, character by character */
    s = (char *) s1;
    s1e = s + (int) ls1 - ls2 + 1;
    while (s < s1e) {
        
        /* Search for first character in pattern string */
        if (*s == cfirst || *s == ocfirst) {
            
            /* If single character search, return */
            if (ls2 == 1)
                return (s);
            
            /* Search for last character in pattern string if first found */
            sl = s[ls2-1];
            if (sl == clast || sl == oclast) {
                
                /* If two-character search, return */
                if (ls2 == 2)
                    return (s);
                
                /* If 3 or more characters, check for rest of search string */
                i = 1;
                while (i < ls2 && (s[i] == (char) s2[i] || s[i] == os2[i]))
                    i++;
                
                /* If entire string matches, return */
                if (i >= ls2) {
                    free (os2);
                    return (s);
                }
            }
	    }
        s++;
	}
    if (os2 != NULL)
        free (os2);
        return (NULL);
}


int
notnum (string)

const char *string;	/* Character string */
{
    if (isnum (string))
        return (0);
    else
        return (1);
}


/* ISNUM-- Return 1 if string is an integer number,
 2 if floating point,
 3 if sexigesimal, with or without decimal point
 else 0
 */

int
isnum (string)

const char *string;	/* Character string */
{
    int lstr, i, nd, cl;
    char cstr, cstr1;
    int fpcode;
    
    /* Return 0 if string is NULL */
    if (string == NULL)
        return (0);
    
    lstr = strlen (string);
    nd = 0;
    cl = 0;
    fpcode = 1;
    
    /* Return 0 if string starts with a D or E */
    cstr = string[0];
    if (cstr == 'D' || cstr == 'd' ||
        cstr == 'E' || cstr == 'e') {
        return (0);
	}
    
    /* Remove trailing spaces */
    while (string[lstr-1] == ' ')
        lstr--;
    
    /* Numeric strings contain 0123456789-+ and d or e for exponents */
    for (i = 0; i < lstr; i++) {
        cstr = string[i];
        if (cstr == '\n')
            break;
        
        /* Ignore leading spaces */
        if (cstr == ' ' && nd == 0)
            continue;
        
        if ((cstr < 48 || cstr > 57) &&
            cstr != '+' && cstr != '-' &&
            cstr != 'D' && cstr != 'd' &&
            cstr != 'E' && cstr != 'e' &&
            cstr != ':' && cstr != '.')
            return (0);
        else if (cstr == '+' || cstr == '-') {
            if (string[i+1] == '-' || string[i+1] == '+')
                return (0);
            else if (i > 0) {
                cstr1 = string[i-1];
                if (cstr1 != 'D' && cstr1 != 'd' &&
                    cstr1 != 'E' && cstr1 != 'e' &&
                    cstr1 != ':' && cstr1 != ' ')
                    return (0);
            }
	    }
        else if (cstr >= 47 && cstr <= 57)
            nd++;
        
        /* Check for colon */
        else if (cstr == 58)
            cl++;
        if (cstr=='.' || cstr=='d' || cstr=='e' || cstr=='d' || cstr=='e')
            fpcode = 2;
	}
    if (nd > 0) {
        if (cl)
            fpcode = 3;
        return (fpcode);
	}
    else
        return (0);
}


/* NUMDEC -- Return number of decimal places in numeric string (-1 if not number) */

int
numdec (string)

const char *string;	/* Numeric string */
{
    char *cdot;
    int lstr;
    
    if (notnum (string) && !strchr (string, ':'))
        return (-1);
    else {
        lstr = strlen (string);
        if ((cdot = strchr (string, '.')) == NULL)
            return (0);
        else
            return (lstr - (cdot - string) - 1);
    }
}


#ifdef USE_SAOLIB
int set_saolib(hstring)
void *hstring;
{
    if( *((int *)hstring) == 142857 )
        use_saolib = 1;
        else
            use_saolib = 0;
            }

#endif


/* Remove exponent, leading #, and/or trailing zeroes, if reasonable */
void
strfix (string, fillblank, dropzero)

char	*string;	/* String to modify */
int	fillblank;	/* If nonzero, fill blanks with underscores */
int	dropzero;	/* If nonzero, drop trailing zeroes */
{
    char *sdot, *s, *strend, *str, ctemp, *slast;
    int ndek, lstr, i;
    
    /* If number, ignore leading # and remove trailing non-numeric character */
    if (string[0] == '#') {
        strend = string + strlen (string);
        str = string + 1;
        strend = str + strlen (str) - 1;
        ctemp = *strend;
        if (!isnum (strend))
            *strend = (char) 0;
        if (isnum (str)) {
            strend = string + strlen (string);
            for (str = string; str < strend; str++)
                *str = *(str + 1);
	    }
        else
            *strend = ctemp;
	}
    
    /* Remove positive exponent if there are enough digits given */
    if (isnum (string) > 1 && strsrch (string, "E+") != NULL) {
        lstr = strlen (string);
        ndek = (int) (string[lstr-1] - 48);
        ndek = ndek + (10 * ((int) (string[lstr-2] - 48)));
        if (ndek < lstr - 7) {
            lstr = lstr - 4;
            string[lstr] = (char) 0;
            string[lstr+1] = (char) 0;
            string[lstr+2] = (char) 0;
            string[lstr+3] = (char) 0;
            sdot = strchr (string, '.');
            if (ndek > 0 && sdot != NULL) {
                for (i = 1; i <= ndek; i++) {
                    *sdot = *(sdot+1);
                    sdot++;
                    *sdot = '.';
                }
            }
	    }
	}
    
    /* Remove trailing zeroes if they are not significant */
    if (dropzero) {
        if (isnum (string) > 1 && strchr (string, '.') != NULL &&
            strsrch (string, "E-") == NULL &&
            strsrch (string, "E+") == NULL &&
            strsrch (string, "e-") == NULL &&
            strsrch (string, "e+") == NULL) {
            lstr = strlen (string);
            s = string + lstr - 1;
            while (*s == '0' && lstr > 1) {
                if (*(s - 1) != '.') {
                    *s = (char) 0;
                    lstr --;
                }
                s--;
            }
	    }
	}
    
    /* Remove trailing decimal point */
    lstr = strlen (string);
    s = string + lstr - 1;
    if (*s == '.')
        *s = (char) 0;
        
    /* Replace embedded blanks with underscores, if requested to */
        if (fillblank) {
            lstr = strlen (string);
            slast = string + lstr;
            for (s = string; s < slast; s++) {
                if (*s == ' ') *s = '_';
            }
        }
    
    return;
    
}

/* Oct 28 1994	New program
 *
 * Mar  1 1995	Search for / after second quote, not first one
 * May  2 1995	Initialize line in HGETC; deal with logicals in HGETL better
 * May  4 1995	Declare STRSRCH in KSEARCH
 * Aug  7 1995  Fix line initialization in HGETC
 * Dec 22 1995	Add HGETRA and HGETDEC to get degrees from xx:xx:xx.xxx string
 *
 * Jan 26 1996	Fix HGETL to not crash when parameter is not present
 * Feb  1 1996	Fix HGETC to deal with quotes correctly
 * Feb  1 1996	Fix HGETDEG to deal with sign correctly
 * Feb  6 1996	Add HGETS to update character strings
 * Feb  8 1996	Fix STRSRCH to find final characters in string
 * Feb 23 1996	Add string to degree conversions
 * Apr 26 1996	Add HGETDATE to get fractional year from date string
 * May 22 1996	Fix documentation; return double from STR2RA and STR2DEC
 * May 28 1996	Fix string translation of RA and Dec when no seconds
 * Jun 10 1996	Remove unused variables after running lint
 * Jun 17 1996	Fix bug which failed to return single character strings
 * Jul  1 1996	Skip sign when reading declination after testing for it
 * Jul 19 1996	Do not divide by 15 if RA header value is already in degrees
 * Aug  5 1996	Add STRNSRCH to search strings which are not null-terminated
 * Aug  6 1996	Make minor changes after lint
 * Aug  8 1996	Fix ksearch bug which finds wrong keywords
 * Aug 13 1996	Fix sign bug in STR2DEC for degrees
 * Aug 26 1996	Drop unused variables ICOL0, NLINE, PREVCHAR from KSEARCH
 * Sep 10 1996	Fix header length setting code
 * Oct 15 1996	Clean up loops and fix ICOL assignment
 * Nov 13 1996	Handle integer degrees correctly in STR2DEC
 * Nov 21 1996	Make changes for Linux thanks to Sidik Isani
 * Dec 12 1996	Add ISNUM to check to see whether strings are numbers
 *
 * Jan 22 1997	Add ifdefs for Eric Mandel (SAOtng)
 * Jan 27 1997	Convert to integer through ATOF so exponents are recognized
 * Jul 25 1997	Implement FITS version of ISO date format
 *
 * Feb 24 1998	Implement code to return IRAF multiple-keyword strings
 * Mar 12 1998	Add subroutine NOTNUM
 * Mar 27 1998	Add changes to match SKYCAT version
 * Apr 30 1998	Add BLSEARCH() to find blank lines before END
 * May 27 1998	Add HGETNDEC() to get number of decimal places in entry
 * Jun  1 1998	Add VMS patch from Harry Payne at StSci
 * Jun 18 1998	Fix code which extracts tokens from string values
 * Jul 21 1998	Drop minus sign for values of -0
 * Sep 29 1998	Treat hyphen-separated date as old format if 2-digit year
 * Oct  7 1998	Clean up search for last blank line
 *
 * Apr  5 1999	Check lengths of strings before copying them
 * May  5 1999	values.h -> POSIX limits.h: MAXINT->INT_MAX, MAXSHORT->SHRT_MAX
 * Jul 15 1999	Add hgetm() options of 1- or 2-digit keyword extensions
 * Oct  6 1999	Add gethlength() to return header length
 * Oct 14 1999	In ksearch(), search only to null not to end of buffer
 * Oct 15 1999	Return 1 from hgetndec() if successful
 * Oct 20 1999	Drop unused variable after lint (val in hgetndec)
 * Dec  3 1999	Fix isnum() to reject strings starting with a d or e
 * Dec 20 1999	Update hgetdate() to get minutes and seconds right
 *
 * Feb 10 2000	Parse RA and Dec with spaces as well as colons as separators
 * Feb 11 2000	Add null at end of multi-line keyword value character string
 * Feb 25 2000	Change max search string length from 57600 to 256000
 * Mar 15 2000	Deal with missing second quotes in string values
 * Mar 17 2000	Return 2 from isnum() if number is floating point (.de)
 * Mar 17 2000	Ignore leading # for numeric values in header
 * Mar 21 2000	Implement -n to get string value starting with nth token
 * Apr  5 2000	Reject +- in isnum()
 * Jun  9 2000	Read keyword values even if no equal sign is present
 * Sep 20 2000	Ignore linefeed at end of number in isnum()
 * Oct 23 2000	Fix handling of embedded + or - in isnum()
 *
 * Jan 19 2000	Return 0 from isnum(), str2ra(), and str2dec() if string is null
 * Mar 30 2001	Fix header length finding algorithm in ksearch()
 * Jul 13 2001	Make val[] static int instead of int; drop unused variables
 * Sep 12 2001	Read yyyy/mm/dd dates as well as dd/mm/yyyy
 * Sep 20 2001	Ignore leading spaces in str2dec()
 * Sep 20 2001	Ignore trailing spaces in isnum()
 *
 * Apr  3 2002	Add hgetr8c(), hgeti4c(), and hgetsc() for multiple WCS handling
 * Apr 26 2002	Fix bug in hgetsc(), hgeti4c(), and hgetr8c() found by Bill Joye
 * Jun 26 2002	Do not drop leading or trailing spaces in multi-line values
 * Aug  6 2002	Add strcsrch() and strncsrch() for case-insensitive searches
 * Aug 30 2002	Fix bug so strcsrch() really is case-insensitive
 * Oct 20 2003	Add numdec() to return number of decimal places in a string
 * Dec  9 2003	Fix numdec() to return 0 if no digits after decimal point
 *
 * Feb 26 2004	Extract value from keyword=value strings within a keyword value
 * Apr  9 2004	Use strncsrch() in ksearch() to find differently-cased keywords
 * Apr 28 2004	Free os2 in strncsrch() only if it is allocated
 * Jul 13 2004	Accept D, d, E, or e as exponent delimiter in floating points
 * Aug 30 2004	Change numdec() to accept sexigesimal numbers (:'s)
 *
 * Jun 27 2005	Drop unused variables
 * Aug 30 2005	Adjust code in hlength()
 *
 * Jun 20 2006	Initialize uninitialized variables in strnsrch()
 * Jun 29 2006	Add new subroutine strfix() to clean strings for other uses
 * Jul 13 2006	Increase maximum number of multiline keywords from 20 to 500
 *
 * Jan  4 2007  Declare header, keyword to be const
 * Jan  4 2007	Change WCS letter from char to char*
 * Feb 28 2007	If header length is not set in hlength, set it to 0
 * May 31 2007	Add return value of 3 to isnum() if string has colon(s)
 * Aug 22 2007	If closing quote not found, make one up
 */
