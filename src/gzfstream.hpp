/*
 * gzfstream.hpp
 *
 *  Created on: Sep 30, 2018
 *      Author: yuankai
 */

#ifndef GZFSTREAM_HPP_
#define GZFSTREAM_HPP_

# include <iosfwd>
# include <boost/iostreams/categories.hpp>
# include <boost/iostreams/positioning.hpp>
# include <boost/iostreams/stream.hpp>
# include <boost/iostreams/concepts.hpp>

# include <zlib.h>

class igzfd : public boost::iostreams::source
{
public:
	typedef char char_type;
	typedef boost::iostreams::source_tag category;

	igzfd(const std::string& path);
	igzfd();

	std::streamsize read(char* s, std::streamsize n);
	boost::iostreams::stream_offset seek(boost::iostreams::stream_offset off,\
			std::ios_base::seekdir way);
	int open(const std::string& path);
	void close();

private:
	gzFile fd;
};

class ogzfd : public boost::iostreams::sink
{
public:
	typedef char char_type;
	typedef boost::iostreams::sink_tag category;

	ogzfd(const std::string& path);
	ogzfd();

	std::streamsize write(const char* s, std::streamsize n);
	boost::iostreams::stream_offset seek(boost::iostreams::stream_offset off,\
			std::ios_base::seekdir way);
	int open(const std::string& path);
	void close();

private:
	gzFile fd;
};

class gzifstream : public boost::iostreams::stream<igzfd>
{
public:
	gzifstream() : boost::iostreams::stream<igzfd>("") {};
	gzifstream(const std::string& path) : boost::iostreams::stream<igzfd>(path) {};
	~gzifstream();

	int open(const std::string& path);
};

class gzofstream : public boost::iostreams::stream<ogzfd>
{
public:
	gzofstream() : boost::iostreams::stream<ogzfd>("") {};
	gzofstream(const std::string& path) : boost::iostreams::stream<ogzfd>(path) {};
	~gzofstream();

	int open(const std::string& path);
};

#endif /* GZFSTREAM_HPP_ */
