/*
 * Timings.hpp
 *
 *  Created on: Dec 1, 2015
 *      Author: fabio
 */

#ifndef TIMINGS_HPP_
#define TIMINGS_HPP_

#include <vector>
#include <numeric>
#include <chrono>
#include <ctime>
#include <cstdarg>

#include <sys/time.h>

namespace tmn {

enum { ST=0, SP=1, GT=2 };

/**
 * Calculate exact elapsed time and return time in msec
 */
static inline double elapsedTime(int tag) {
	static struct timeval tv_start = {0,0};
	static struct timeval tv_stop  = {0,0};

	double res=0.0;
	switch(tag) {
	case ST:{
		gettimeofday(&tv_start,NULL);
	} break;
	case SP:{
		gettimeofday(&tv_stop,NULL);
		long sec  = (tv_stop.tv_sec  - tv_start.tv_sec);
		long usec = (tv_stop.tv_usec - tv_start.tv_usec);

		if(usec < 0) {
			usec += 1000000;
		}
		res = ((double)(sec*1000)+ ((double)usec)/1000.0);
	} break;
	case GT: {
		long sec  = (tv_start.tv_sec  - tv_stop.tv_sec);
		long usec = (tv_start.tv_usec - tv_stop.tv_usec);

		if(usec < 0) {
			--sec;
			usec += 1000000;
		}
		res = ((double)(sec*1000)+ ((double)usec)/1000.0);
	} break;
	default:
		res=0;
		break;
	}
	return res;
}

inline std::string timeStamp() {
	char buf[50];
	char rm[] = ":-";

	std::time_t t = time(0);
	struct tm * now = localtime( & t );
	strftime(buf, sizeof(buf), "%Fh%T", now);

	std::string tmstmp(buf);
	for(int i=0; i<2; ++i)
		tmstmp.erase( std::remove(tmstmp.begin(), tmstmp.end(), rm[i]),
				tmstmp.end() );

	return tmstmp;
}

inline std::string timeStamp2() {
	char buf[50];

	std::time_t t = time(0);
	struct tm * now = localtime( & t );
	strftime(buf, sizeof(buf), "%F", now);

	std::string tmstmp(buf);
	//std::replace(tmstmp.begin(), tmstmp.end(), ':', '-');

	return tmstmp;
}

} // namespace tmn


#define fw(what) std::forward<decltype(what)>(what)

/**
 * @ class measure
 * @ brief Class to measure the execution time of a callable
 *
 * Taken from Nikos Athanasiou
 * (https://github.com/picanumber/bureaucrat)
 */
template <typename TimeT = std::chrono::milliseconds,
		class ClockT = std::chrono::system_clock >
struct Measure {

	/**
	 * @ fn    execution
	 * @ brief Returns the quantity (count) of the elapsed time as TimeT units
	 */
	template<typename F, typename ...Args>
	static typename TimeT::rep execution(F&& func, Args&&... args) {
		auto start = ClockT::now();
		fw(func)(std::forward<Args>(args)...);
		auto duration = std::chrono::duration_cast<TimeT>(ClockT::now() - start);

		return duration.count();
	}

	/**
	 * @ fn    duration
	 * @ brief Returns the duration (in chrono's type system) of the elapsed time
	 */
	template<typename F, typename... Args>
	static TimeT duration(F&& func, Args&&... args)	{
		auto start = ClockT::now();
		fw(func)(std::forward<Args>(args)...);

		return std::chrono::duration_cast<TimeT>(ClockT::now() - start);
	}
};

/*
 * Measure - Sample usage
 *
int main() {
	std::vector<int> ar(2000000);
	std::iota(begin(ar), end(ar), 0); // fills *ar from 0 to *ar.size()-1

	// 1. Client just wants a timing result : call execution
	std::cout << measure<std::chrono::nanoseconds>::execution(
			std::accumulate<decltype(begin(ar)), int>, begin(ar), end(ar), 0) << std::endl;

	// 2. Client want to preprocess timing : call duration, process result and then query the count
	auto avg = measure<>::duration(std::accumulate<decltype(begin(ar)), int>, begin(ar), end(ar), 0);
	avg     += measure<>::duration(std::accumulate<decltype(begin(ar)), int>, begin(ar), end(ar), 0);
	avg     += measure<>::duration(std::accumulate<decltype(begin(ar)), int>, begin(ar), end(ar), 0);
	std::cout << (avg / 3.).count() << std::endl;
}
*/



#endif /* TIMINGS_HPP_ */
