#include "handle_io.h"

std::string select_linemark(msgcode c)
{
	std::string linemark;
	switch (c)
	{
	case MSG_EMPTY:
		break;
	case MSG_INFO:
		linemark = "[Info] ";
		break;
	case MSG_WARNING:
		linemark = "[Warning] ";
		break;
	case MSG_RESULT:
		linemark = "[Result] ";
		break;
	case MSG_ERROR:
		linemark = "[Error] ";
		break;
	default:
		;
	}
	return linemark;
}

void handle_message(msgcode c, std::string msg)
{
	if (c == MSG_ERROR)
	{
		std::cin.get();
		std::cout << select_linemark(c) << msg << std::endl;
		std::cout << "Program interrupted. Press any key to exit...";
		std::cin.get();
		exit(1);
	}
	else
		std::cout << select_linemark(c) << msg << std::endl;
}

void handle_message(msgcode c, boost::format msg)
{
	if (c == MSG_ERROR)
	{
		std::cin.get();
		std::cout << select_linemark(c) << msg << std::endl;
		std::cout << "Program interrupted. Press any key to exit...";
		std::cin.get();
		exit(1);
	}
	else
		std::cout << select_linemark(c) << msg << std::endl;
}
