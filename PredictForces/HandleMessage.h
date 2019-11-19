#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <boost/format.hpp>

#include <Eigen/Dense>

enum msgcode : uint8_t {
	MSG_EMPTY,
	MSG_INFO,
	MSG_WARNING,
	MSG_RESULT,
	MSG_ERROR
};

const char interruptInfo[47] = "Program interrupted. Press any key to exit...";

// Matrix formats
const Eigen::IOFormat CleanFmt = Eigen::IOFormat(4, 0, ", ", "\n", "[", "]");

class HandleMessage
{
public:
	HandleMessage();
	HandleMessage(std::string logPath);
	~HandleMessage();

	void handle_message(msgcode c, std::string msg);
	void handle_message(msgcode c, boost::format msg);

	static void print(msgcode c, std::string msg, std::ostream & out = std::cout);
	static void print(msgcode c, boost::format msg, std::ostream & out = std::cout);
	static void print_matrix(Eigen::MatrixXd mat, Eigen::IOFormat format = CleanFmt, std::ostream & out = std::cout);

private:
	static std::string select_linemark(msgcode c);

	std::ofstream targetStream;
};

