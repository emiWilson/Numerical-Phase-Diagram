#ifndef STRINGCONVERT_H
#define STRINGCONVERT_H

#include <string>
#include <sstream>
#include <iomanip>  

bool startsWith(const std::string str, const std::string compare) {
	if(str.length() < compare.length())
		return false;
	else 
		return str.compare(0, compare.length(), compare) == 0;
}

std::string trim(const std::string& str, const std::string& whitespace = " \t")
{
    const unsigned int strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const int strEnd = str.find_last_not_of(whitespace);
    const int strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

template <class T> std::string to_string (const T& t)
{
	std::stringstream ss;
	ss << t;
	return ss.str();
}

template <class T> std::string to_string (const T& t, int width)
{
	std::stringstream ss;
	ss << std::setfill ('0') << std::setw (width);
	ss << t;
	return ss.str();
}

template <class T> T string_convert(const std::string& s) {
	std::istringstream i(s);
	T x;
	if (!(i >> x)) {
		//return NAN;
	}
   return x;
}



#endif