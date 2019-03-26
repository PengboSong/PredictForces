#pragma once
#include <algorithm> 
#include <cctype>
#include <locale>
#include <string>
#include <sstream>

static void ltrim(std::string &s);

static void rtrim(std::string &s);

static void trim(std::string &s);

std::string sslice(size_t begin, size_t end, std::string &in);

int ssliceti(size_t begin, size_t end, std::string &in);

double sslicetd(size_t begin, size_t end, std::string &in);
