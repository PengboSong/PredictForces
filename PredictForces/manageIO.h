#pragma once
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void handle_info(string info);

void handle_warning(string info);

void handle_error(string info);

void handle_result(string info);

void handle_result(string info, vector<string> buf);
