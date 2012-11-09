//
// Copyright (c) 2009-2010  Shuntaro Yamazaki (shun-yamazaki (at) aist.go.jp)
//  and the National Institute of Advanced Industrial Science and Technology
//
// $Id: MiscUtil.h 5727 2012-03-12 17:48:51Z shun $
//

#pragma once

#ifndef TRACE
#include <cstdio>
#define TRACE printf
#endif

#ifndef _ASSERTE
#include <cassert>
#define _ASSERTE assert
#endif

#include <sys/types.h>
#include <sys/stat.h>

#undef min
#undef max

namespace slib
{

inline
void ThrowRuntimeError(const char *fmt, ...)
{
	va_list param;
	va_start(param, fmt);
	char message[1024];
	_vsnprintf(message, 1023, fmt, param);
	message[1023] = 0;
	TRACE(message);
	TRACE("\n");
	throw std::runtime_error(message);
}

inline
void ThrowLogicError(const char *fmt, ...)
{
	va_list param;
	va_start(param, fmt);
	char message[1024];
	_vsnprintf(message, 1023, fmt, param);
	message[1023] = 0;
	TRACE(message);
	throw std::logic_error(message);
}

inline
std::string string_format(const char *fmt, ...)
{
	std::string buffer;
	va_list     arguments;

	_ASSERTE(fmt);
	va_start(arguments, fmt);
	try
	{
		int length = _vscprintf(fmt, arguments);
		if (length < 0)
		{
			throw std::runtime_error("_vscprintf() failed.");
		}
		buffer.assign(length, 0);
		int result = vsprintf(&buffer[0], fmt, arguments);
		if (result < 0)
		{
			throw std::runtime_error("vsprintf() failed.");
		}
	}
	catch (...)
	{
		va_end(arguments);
		throw;
	}
	va_end(arguments);
	return buffer;
}

inline
const std::string string_replace(const std::string& str, const std::string& oldstr, const std::string& newstr)
{
	if (str.empty() || oldstr.empty())
	{
		return str;
	}
	std::string result = str;
	for (std::string::size_type pos(0); (pos = result.find(oldstr, pos)) != std::string::npos; pos += newstr.length())
	{
		result.replace(pos, oldstr.length(), newstr);
	}
	return result;
}

inline
const std::string string_trim(const std::string& str, const char *delim = " \n\r\t")
{
	const int p1(str.find_first_not_of(delim));
	if (p1 == std::string::npos)
	{
		return std::string();
	}
	const int p2(str.find_last_not_of(delim));
	return str.substr(p1, p2 - p1 + 1);
}

/*
inline
const std::vector<std::string> string_token(const std::string& str)
{
	std::istringstream is(str);
  std::vector<std::string> tokens;
  std::string token;
  for(;;){
    token.clear();
    is >> token;
    if(token.empty()){
      break;
    }
    tokens.push_back(token);
  }
}
*/
inline
bool IsPowerOfTwo(const int n)
{
	return !(n & (n - 1));
}

inline
int GetLeastPowerOfTwo(const int n)
{
	int num = n;
	while (!IsPowerOfTwo(num))
	{
		num++;
	}
	return num;
}

inline
bool file_exists(const std::string& filename)
{
	struct _stat	 buf;
	return !_stat(filename.c_str(), &buf);
}

} // slib
