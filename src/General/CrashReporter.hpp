#ifndef CRASH_REPORTER_HPP
#define CRASH_REPORTER_HPP

#include <string>
#include <fstream>
#include <iostream>

enum CRASH_TYPE {
	CT_TEXTURE = 1,
	CT_SOUND,
	CT_BORDER,
	CT_TRIGGER,
	CT_NEED_KEY
};

class CrashReporter : public std::exception {
    public:
        CrashReporter (const std::string &inf);
        ~CrashReporter ();
        void Report (const std::string &file);
		const char* what() const noexcept;
    protected:
        std::string inf;
        //const std::string reportFile = "2.txt";
};

CrashReporter::CrashReporter (const std::string &inf) : inf(inf) {}

CrashReporter::~CrashReporter () {}

void CrashReporter::Report (const std::string &reportfile) {
	if (reportfile == "") {
 		std::cerr << inf;
 	} else {
		std::string message;
		time_t date = time(0);
		message += ctime(&date);
		message.insert(message.end() - 1, ':');

 		std::ofstream file;
 		file.open(reportfile);
 		file << message << "\n\n" << inf << "\n";
 		file.close();
 	}
}

// void CrashReporter::Report (CRASH_TYPE type, const std::string &path) {
// 	std::string message;
// 	time_t date = time(0);
// 	message += ctime(&date);
// 	message.insert(message.end() - 1, ':');
// 	switch (type) {
// 		case CT_TEXTURE:
// 			message += "Texture error: can't open texture in \"" + path + "\" in object " + name;
// 			break;
// 		case CT_SOUND:
// 			message += "Sound error: can't load sound in \"" + path + "\" in object " + name;
// 			break;
// 		case CT_BORDER:
// 			message += "Border error: can't";
// 			break;
// 		case CT_TRIGGER:
// 			message += "Trigger error: can't";
// 			break;
// 		default:
// 			message += "Unknown error: not registred error in object " + name;
// 			break;
// 	}
// 	message += "\n";
// 	if (reportFile == "") {
// 		std::cerr << message;
// 	} else {
// 		std::ofstream file;
// 		file.open(reportFile);
// 		file << message;
// 		file.close();
// 	}
// 	exit(-type);
// }

const char* CrashReporter::what() const noexcept {
	return inf.c_str();
}

#endif