/*
	This module contains misc routines.

	Kalev Kask, March 2002.
*/

#include <stdlib.h>
#include "Utils/MiscUtils.hxx"

#if defined WINDOWS || _WINDOWS
#include <windows.h>
#endif

#include <cstring>
#include <string>
#include <time.h>
#include <errno.h>

#ifdef LINUX
#include <chrono>
#endif


INT64 ARE::GetTimeInMilliseconds(void)
{
#if defined WINDOWS || _WINDOWS
	FILETIME ft ;
	GetSystemTimeAsFileTime(&ft) ;
	__int64 t = ft.dwHighDateTime ;
	t <<= 32 ;
	t += ft.dwLowDateTime ;
	// 64-bit value representing the number of 100-nanosecond intervals since January 1, 1601 (UTC).
	// convert to milliseconds
	t /= 10000 ;
	return t ;
#else
	return std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
#endif
}


int ARE::ExtractVarValuePairs(char *BUF, int L, std::list<std::pair<std::string,std::string> > & AssignmentList)
{
	AssignmentList.clear() ;
	if (NULL == BUF)
		return 0 ;
	if (L < 0)
		L = strlen(BUF) ;

	int i, j, n = 0 ;
	char *s = BUF ;
	std::string sN, sV ;
	for (i = 0 ; i < L ; i++) {
		if ('\n' != BUF[i]) continue ;
		int l = BUF + i - s ;
		if (l > 0) {
			// [s, s+l) is a string
			for (j = 0 ; j < l ; j++)
				{ if ('=' == s[j]) break ; }
			if (j < l) {
				const char *sName = s ;
				const char *sValue = s + j + 1 ;
				int lN = j ;
				int lV = l - j - 1 ;
				if (lN > 0) sN.assign(sName, lN) ; else sN.erase() ;
				if (lV > 0) sV.assign(sValue, lV) ; else sV.erase() ;
				std::pair<std::string,std::string> assignment(sN,sV) ;
				AssignmentList.push_back(assignment) ;
				}
			}
		if ('\r' == BUF[i+1])
			i++ ;
		s = BUF + i + 1 ;
		}
	int l = BUF + i - s ;
	if (l > 0) {
		// [s, s+l) is a string
		for (j = 0 ; j < l ; j++)
			{ if ('=' == s[j]) break ; }
		if (j < l) {
			const char *sName = s ;
			const char *sValue = s + j + 1 ;
			int lN = j ;
			int lV = l - j - 1 ;
			if (lN > 0) sN.assign(sName, lN) ; else sN.erase() ;
			if (lV > 0) sV.assign(sValue, lV) ; else sV.erase() ;
			std::pair<std::string,std::string> assignment(sN,sV) ;
			AssignmentList.push_back(assignment) ;
			}
		}

	return n ;
}


int ARE::ExtractParameterValue(/* IN */ std::string & Paramater, std::list<std::pair<std::string,std::string> > AssignmentList, /* OUT */ std::string & Value)
{
	Value.erase() ;
	for (std::list<std::pair<std::string,std::string> >::iterator i = AssignmentList.begin() ; i != AssignmentList.end() ; i++) {
		std::pair<std::string,std::string> & assignment = *i ;
		std::string & sN = assignment.first ;
		if (Paramater == sN)
			{ Value = assignment.second ; return 0 ; }
		}

	return 1 ;
}


int ARE::ExtractTokens(const std::string & S, std::list<std::string> & Tokens, int bIncludeEmptyTokens, char wcSeparator)
{
//	std::list<std::wstring>::iterator it1 = Tokens.begin(), it2 ;
	Tokens.erase(Tokens.begin(), Tokens.end()) ;

	// RULES:
	// * ';' is a default separator between tokens; but a different separator symbol could be given (e.g. ',').
	// * '{', '}', '(', ')' allow nested tokens.
	// * anything between double quotation marks ('"') is a string; these double quotation marks are not part of value of the string.
	//   i.e. at some point when the token is processed, these double quotation marks will be stripped away.
	// * in order to include a double quotation mark in a string, it should be doubled; i.e. '""' is equals '"' such that
	//   it is considered part of the string (i.e. not beginning/end of a string).

	int embedCount = 0 ;
	int pS = 0 ;
	bool insideString = false, insideBlindString = false ;
	std::string sEmpty ;
	char delimiter = 0 ;
	std::string s ;
	unsigned int i ;
	for (i = 0 ; i < S.length() ; i++) {
		const wchar_t & tc = S[i] ;

		// if currently inside a string, process here.
		if (insideString) {
			if ('"' == tc) {
				int ip1 = i+1 ;
				if (S[ip1]) {
					if ('"' == S[ip1]) {
						++i ;
						continue ;
						}
					}
				insideString = false ;
				}
			continue ;
			}
		if (insideBlindString) {
			if (delimiter == tc)
				insideBlindString = false ;
			continue ;
			}

		// check if entering a string
		if ('"' == tc) {
			insideString = true ;
			continue ;
			}
		if ('@' == tc) {
			// we have "@<delimiter>...<delimiter>"; there must be at least 2 more characters
			delimiter = S[++i] ;
			insideBlindString = true ;
			continue ;
			}

		// process rest.
		if ('}' == tc || ')' == tc) {
			--embedCount ;
			continue ;
			}
		else if ('{' == tc || '(' == tc) {
			++embedCount ;
			continue ;
			}
		// skip everything inside '{','}', '(', ')'
		if (embedCount > 0)
			continue ;
		if (wcSeparator == tc) {
			int l = i - pS ;
			if (l > 0) {
				// 2008-07-14 KK : while tokens may be escaped/double-quoted, don't un-escape/un-double-quote them here.
				// we may not know when/how to do it correctly.
				// e.g. complete enum-list should not be done here : "OK"{"OK";"Problem"}
				Tokens.push_back(S.substr(pS, l)) ;
/*
				if (0 == __AllocateString(__MSstring_pStemp, __MSstring_StempL, l+1)) {
//					s = S.substr(pS, l) ;
					int OUTlen = 0 ;
					CleanupString(S.c_str()+pS, l, __MSstring_pStemp, OUTlen) ;
					if (OUTlen > 0) {
						s.assign(__MSstring_pStemp, OUTlen) ;
						Tokens.Insert(s) ;
						}
					else if (bIncludeEmptyTokens)
						Tokens.Insert(sEmpty) ;
					}
*/
				}
			else if (bIncludeEmptyTokens)
				Tokens.push_back(sEmpty) ;
			pS = i+1 ;
			}
		}
	// check for errors
	if (insideString || insideBlindString) {
		// error: S ends before the nested string ('"' or special delimiter)
		Tokens.erase(Tokens.begin(), Tokens.end()) ;
		return 1 ;
		}
	int l = i - pS ;
	if (l > 0) {
		// 2008-07-14 KK : while tokens may be escaped/double-quoted, don't un-escape/un-double-quote them here.
		// we may not know when/how to do it correctly.
		// e.g. complete enum-list should not be done here : "OK"{"OK";"Problem"}
		Tokens.push_back(S.substr(pS, l)) ;
/*
		if (0 == __AllocateString(__MSstring_pStemp, __MSstring_StempL, l+1)) {
			int OUTlen = 0 ;
			CleanupString(S.c_str()+pS, l, __MSstring_pStemp, OUTlen) ;
			if (OUTlen > 0) {
				s.assign(__MSstring_pStemp, OUTlen) ;
				Tokens.Insert(s) ;
				}
			else if (bIncludeEmptyTokens)
				Tokens.Insert(sEmpty) ;
			}
*/
		}
	else if (bIncludeEmptyTokens)
		Tokens.push_back(sEmpty) ;

	return 0 ;
}
