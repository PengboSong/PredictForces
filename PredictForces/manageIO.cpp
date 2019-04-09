#include "manageIO.h"

void handle_info(string info)
{
	cout << "[Info] " << info << endl;
}

void handle_warning(string info)
{
	cout << "[Warning] " << info << endl;
}

void handle_error(string info)
{
	cout << "[Error] " << info << endl;
	cout << "Program interrupted. Press any key to exit...";
	cin.get();
}

void handle_result(string info)
{
	cout << "[Result] " << info << endl;
}

void handle_result(string info, vector<string> buf)
{
	cout << "[Result] " << info << endl;
	for (vector<string>::iterator it = buf.begin(); it != buf.end(); ++it)
		cout << *it << endl;
}
