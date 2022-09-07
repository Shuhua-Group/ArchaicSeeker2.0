/*
 * gzfstream.cpp
 *
 *  Created on: Sep 30, 2018
 *      Author: yuankai
 */

# include "gzfstream.hpp"

igzfd::igzfd(const std::string& path)
{
	fd = gzopen(path.c_str(), "rb");
}

std::streamsize igzfd::read(char* s, std::streamsize n)
{
	std::streamsize rd = gzread(fd, s, n);
	if(rd != 0)
		return rd;
	else
		return -1;
}

boost::iostreams::stream_offset igzfd::seek(boost::iostreams::stream_offset off, \
		std::ios_base::seekdir way)
{
	if(off < 0)
		throw std::ios_base::failure("bad seek offset");
	if(way == std::ios_base::cur)
		gzseek(fd, off, SEEK_CUR);
	else if(way == std::ios_base::beg)
		gzseek(fd, off, SEEK_SET);
	else
		throw std::ios_base::failure("bad seek direction");
	return gztell(fd);
}

int igzfd::open(const std::string& path)
{
	fd = gzopen(path.c_str(), "rb");
	int ret = 1;
	if(NULL == fd)
		ret = 0;
	return ret;
}

void igzfd::close()
{
	gzclose(fd);
}

ogzfd::ogzfd(const std::string& path)
{
	fd = gzopen(path.c_str(), "wb");
}

std::streamsize ogzfd::write(const char* s, std::streamsize n)
{
	std::streamsize wd = gzwrite(fd, s, n);
	return wd;
}

boost::iostreams::stream_offset ogzfd::seek(boost::iostreams::stream_offset off, \
		std::ios_base::seekdir way)
{
	if(off < 0)
		throw std::ios_base::failure("bad seek offset");
	if(way == std::ios_base::cur)
		gzseek(fd, off, SEEK_CUR);
	else if(way == std::ios_base::beg)
		gzseek(fd, off, SEEK_SET);
	else
		throw std::ios_base::failure("bad seek direction");
	return gztell(fd);
}

int ogzfd::open(const std::string& path)
{
	fd = gzopen(path.c_str(), "wb");
	int ret = 1;
	if(NULL == fd)
		ret = 0;
	return ret;
}

void ogzfd::close()
{
	gzclose(fd);
}

gzifstream::~gzifstream()
{
	this->operator->()->close();
}

int gzifstream::open(const std::string& path)
{
	return this->operator->()->open(path);
}

gzofstream::~gzofstream()
{
	this->operator->()->close();
}

int gzofstream::open(const std::string& path)
{
	return this->operator->()->open(path);
}



