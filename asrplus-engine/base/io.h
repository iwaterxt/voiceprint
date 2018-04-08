/*
 * io.h
 *
 *  Created on: Feb 1, 2018
 *      Author: tao
 */

#ifndef IO_H_
#define IO_H_

# include <iostream>
# include <sstream>

void CheckToken(const char *token) ;
void ReadToken(std::istream &is, bool binary, std::string *token);

void ExpectToken(std::istream &is, bool binary, const char *token);
void ExpectToken(std::istream &is, bool binary, const std::string & token);

void ReadBasicType(std::istream &is, bool binary, float *f);
void ReadBasicType(std::istream &is, bool binary, int *f);
void ReadBasicType(std::istream &is, bool binary, bool *b);

#endif /* IO_H_ */
