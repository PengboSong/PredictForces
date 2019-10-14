#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <boost/format.hpp>

enum msgcode : uint8_t {
	MSG_EMPTY,
	MSG_INFO,
	MSG_WARNING,
	MSG_RESULT,
	MSG_ERROR
};

std::string select_linemark(msgcode c);

void handle_message(msgcode c, std::string msg);
void handle_message(msgcode c, boost::format msg);
