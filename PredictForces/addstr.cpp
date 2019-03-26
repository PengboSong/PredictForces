#include "addstr.h"

// trim from start (in place)
static inline void ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
		return !std::isspace(ch);
	}));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
		return !std::isspace(ch);
	}).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
	ltrim(s);
	rtrim(s);
}

std::string sslice(size_t begin, size_t end, std::string &in)
{
	std::string out = "";
	if (begin < in.length() && end < in.length())
	{
		if (end > begin)
			for (size_t i = begin; i < end; i++)
				out += in[i];
		else if (end < begin)
			for (size_t i = begin; i < end; i--)
				out += in[i];
	}
	trim(out);
	return out;
}

int ssliceti(size_t begin, size_t end, std::string &in)
{
	int value = 0;
	std::stringstream tmp;
	if (begin < in.length() && end < in.length())
	{
		if (end > begin)
			for (size_t i = begin; i < end; i++)
				tmp << in[i];
		else if (end < begin)
			for (size_t i = begin; i < end; i--)
				tmp << in[i];
	}
	tmp >> value;
	return value;
}

double sslicetd(size_t begin, size_t end, std::string &in)
{
	double value = 0.0;
	std::stringstream tmp;
	if (begin < in.length() && end < in.length())
	{
		if (end > begin)
			for (size_t i = begin; i < end; i++)
				tmp << in[i];
		else if (end < begin)
			for (size_t i = begin; i < end; i--)
				tmp << in[i];
	}
	tmp >> value;
	return value;
}
