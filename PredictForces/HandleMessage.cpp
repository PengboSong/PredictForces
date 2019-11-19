#include "HandleMessage.h"

HandleMessage::HandleMessage()
{
}

HandleMessage::HandleMessage(std::string logPath)
{
	targetStream.open(logPath, std::ios::app);
}

HandleMessage::~HandleMessage()
{
	if (targetStream.is_open())
		targetStream.close();
}

void HandleMessage::handle_message(msgcode c, std::string msg)
{
	if (c == MSG_ERROR)
	{
		std::cin.get();

		if (targetStream.is_open())
		{
			targetStream << select_linemark(c) << msg << std::endl;
			targetStream.close();
		}

		std::cout << select_linemark(c) << msg << std::endl << interruptInfo;
		std::cin.get();
		exit(1);
	}
	else
	{
		if (targetStream.is_open())
			targetStream << select_linemark(c) << msg << std::endl;
		std::cout << select_linemark(c) << msg << std::endl;
	}
}

void HandleMessage::handle_message(msgcode c, boost::format msg)
{
	if (c == MSG_ERROR)
	{
		std::cin.get();

		if (targetStream.is_open())
		{
			targetStream << select_linemark(c) << msg << std::endl;
			targetStream.close();
		}

		std::cout << select_linemark(c) << msg << std::endl << interruptInfo;
		std::cin.get();
		exit(1);
	}
	else
	{
		if (targetStream.is_open())
			targetStream << select_linemark(c) << msg << std::endl;
		std::cout << select_linemark(c) << msg << std::endl;
	}
}

void HandleMessage::print(msgcode c, std::string msg, std::ostream & out)
{
	if (c == MSG_ERROR)
	{
		std::cin.get();
		out << select_linemark(c) << msg << std::endl << interruptInfo;
		std::cin.get();
		exit(1);
	}
	else
		out << select_linemark(c) << msg << std::endl;
}

void HandleMessage::print(msgcode c, boost::format msg, std::ostream & out)
{
	if (c == MSG_ERROR)
	{
		std::cin.get();
		out << select_linemark(c) << msg << std::endl << interruptInfo;
		std::cin.get();
		exit(1);
	}
	else
		out << select_linemark(c) << msg << std::endl;
}

void HandleMessage::print_matrix(Eigen::MatrixXd mat, Eigen::IOFormat format, std::ostream & out)
{
	out << mat.format(format) << std::endl;
}

std::string HandleMessage::select_linemark(msgcode c)
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
