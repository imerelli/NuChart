/*
 * FileIO.hpp
 *
 *  Created on: Dec 1, 2015
 *      Author: fabio
 */

#ifndef FILEIO_HPP_
#define FILEIO_HPP_

#include "MemoryMapped.h"
#include "Timings.hpp"

#include <string>
#include <cstdio>
#include <cstdlib>

#include <unistd.h>
#include <curl/curl.h>
#include <sys/stat.h>
#include <dirent.h>
#include <fcntl.h>
#include <errno.h>

namespace {

void eraseTrailingWhiteSpaces(std::string& str) {
	std::string whitespaces (" \t\f\v\n\r");

	std::size_t found = str.find_last_not_of(whitespaces);
	if (found!=std::string::npos)
		str.erase(found+1);
	else
		str.clear();            // str is all whitespace
}

// check if file exists
bool exist_f(const char *file) {
#if (defined(_MSC_VER) || defined(__INTEL_COMPILER)) && defined(_WIN32)
	std::fstream fs(file);
	return fs.good();
#else
	struct stat buf;
	return (stat(file, &buf) == 0);
#endif
}

// count lines in a text file (like wc -l does)
inline int countLines(std::ifstream *inf) {
	//std::ifstream inf(path, std::ifstream::in);
	std::string l;
	int n_lines = 0;

	while (std::getline(*inf, l))
		n_lines++;

	return n_lines;
}

// count lines in a text file mapped to memory(like wc -l does)
inline int countLines(MemoryMapped *mmFile) {
	size_t data_sz = mmFile->size();

	long m_numLines = 0;
	auto bff = mmFile->getData();
	auto l = bff + data_sz;
	while (bff && bff!=l) {
		if ((bff = static_cast<const unsigned char*>(memchr(bff, '\n', l-bff))))
			++m_numLines;
		++bff;
	}

	return m_numLines;
}

// count lines in a text file (like wc -l does)
inline int countLines(std::string &inf) {
	int n_lines = 0;
	static const auto BUFFER_SIZE = 256;
	char buf[BUFFER_SIZE+1];
	int fd = ::open(inf.c_str(), O_RDONLY);
	if (fd == -1)
		printf("[ERROR] Open failed\n");

	int advice = POSIX_FADV_SEQUENTIAL | POSIX_FADV_DONTNEED;
	posix_fadvise(fd, 0, 0, advice);  // FDADVICE_SEQUENTIAL

	while(size_t bytes_read = ::read(fd, buf, BUFFER_SIZE)) {
		if(bytes_read == (size_t)-1)
			printf("[ERROR] read failed\n");
		if (!bytes_read)
			break;
		for(char *p = buf; (p = (char*) ::memchr(p, '\n', (buf + bytes_read) - p)); ++p)
			++n_lines;
	}
	::close(fd);

	return n_lines;
}


// writes in vector *files all files contained in the specified directory *dir
// skips hidden files
int getdirFiles(std::string dir, std::vector<std::string> &files) {
    DIR *dp;
    std::string fPath, f;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
    	std::cout << "Error(" << errno << ") opening " << dir << std::endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	fPath = dir;
    	f = std::string(dirp->d_name);

    	if(f.at(0) == '.') continue;

    	fPath.append(f);
        files.push_back(fPath);
        fPath.clear(); f.clear();
    }
    closedir(dp);
    return 0;
}


/******************************************************************************
 * This is an example showing how to get a single file from an FTP server.
 * It delays the actual destination file creation until the first write
 * callback so that it won't create an empty file in case the remote file
 * doesn't exist or something else fails.
 */

struct FtpFile {
	const char *filename;
	FILE *stream;
};

//
static size_t writeFile(void *buffer, size_t size, size_t nmemb, void *stream)
{
	struct FtpFile *out = (struct FtpFile *)stream;
	if(out && !out->stream) {
		/* open file for writing */
		out->stream = fopen(out->filename, "wb");
		if(!out->stream) {
			perror("ERROR, Cannot open file");
			return -1; /* failure, can't open file to write */
		}
	}
	return fwrite(buffer, size, nmemb, out->stream);
}

// download dataset from server and save it locally
// TODO: this is only for BED files. generalizzare
int getGFfile(std::string chr) {
	CURL *curl;
	CURLcode res;

	chr.append(".txt");

	std::string gf_path = "./extdata/human/GF/";
	std::string gf = "gf_";

	gf_path.append(gf).append(chr);

	std::string nu_server = "ftp://fileserver.itb.cnr.it/nuchart/GF/human/";
	nu_server.append(gf).append(chr);

	struct FtpFile ftpfile = { gf_path.c_str(), NULL };

	curl_global_init(CURL_GLOBAL_DEFAULT);
	curl = curl_easy_init();

	if(curl) {
		curl_easy_setopt(curl, CURLOPT_URL, nu_server.c_str());

		// callback to write data
		curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, writeFile);
		curl_easy_setopt(curl, CURLOPT_WRITEDATA, &ftpfile);

		res = curl_easy_perform(curl);

		// cleanup
		curl_easy_cleanup(curl);

		if(CURLE_OK != res) {
			fprintf(stderr, "curl told us %d\n", res);
			return -1;
		} else {
			std::cerr << "File " << chr << " downloaded in \'" << gf_path << "\'" << std::endl;
		}
	} else {
		std::cerr << "Something went wrong. Cannot use curl functions" << std::endl;
		return -1;
	}

	if(ftpfile.stream)
		fclose(ftpfile.stream);

	curl_global_cleanup();

	return 1;
}
// ----------------------------------------------------------------------------

} // namespace



#endif /* FILEIO_HPP_ */
