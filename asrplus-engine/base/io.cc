/*
 * io.cc
 *
 *  Created on: Feb 1, 2018
 *      Author: tao
 */

#include "base/io.h"
#include "base/common.h"
#include <assert.h>
#include <string.h>
#include <limits>
#include <vector>

void CheckToken(const char *token) {
  if (*token == '\0')
    std::cerr << "Token is empty (not a valid token)";
  const char *orig_token = token;
  while (*token != '\0') {
    if (::isspace(*token))
      std::cerr << "Token is not a valid token (contains space): '"
                << orig_token << "'";
    token++;
  }
}

void ReadToken(std::istream &is, bool binary, std::string *str){

	  assert(str != NULL);
	  if (!binary) is >> std::ws;  // consume whitespace.
	  is >> *str;
	  if (is.fail()) {
	    std::cerr << "ReadToken, failed to read token at file position "
	              << is.tellg();
	  }
	  if (!isspace(is.peek())) {
	    std::cerr << "ReadToken, expected space after token, saw instead "
	              << static_cast<char>(is.peek())
	              << ", at file position " << is.tellg();
	  }
	  is.get();  // consume the space.

}

void ExpectToken(std::istream &is, bool binary, const char *token){

	  int pos_at_start = is.tellg();
	  assert(token != NULL);
	  CheckToken(token);  // make sure it's valid (can be read back)
	  if (!binary) is >> std::ws;  // consume whitespace.
	  std::string str;
	  is >> str;
	  is.get();  // consume the space.
	  if (is.fail()) {
	    std::cerr << "Failed to read token [started at file position "
	              << pos_at_start << "], expected " << token;
	  }
	  if (strcmp(str.c_str(), token) != 0) {
	    std::cerr << "Expected token \"" << token << "\", got instead \""
	              << str <<"\".";
	  }
}

void ExpectToken(std::istream &is, bool binary, const std::string & token){

	ExpectToken(is, binary, token.c_str());


}


void ReadBasicType(std::istream &is, bool binary, float *f) {
  assert(f != NULL);
  if (binary) {
    float d;
    int c = is.peek();
    if (c == sizeof(*f)) {
      is.get();
      is.read(reinterpret_cast<char*>(f), sizeof(*f));
    } else if (c == sizeof(d)) {
      ReadBasicType(is, binary, &d);
      *f = d;
    } else {
      std::cerr << "ReadBasicType: expected float, saw " << is.peek()
                << ", at file position " << is.tellg();
    }
  } else {
    is >> *f;
  }
  if (is.fail()) {
    std::cerr << "ReadBasicType: failed to read, at file position "
              << is.tellg();
  }
}

void ReadBasicType(std::istream &is, bool binary, bool *b) {
  assert(b != NULL);
  if (!binary) is >> std::ws;  // eat up whitespace.
  char c = is.peek();
  if (c == 'T') {
      *b = true;
      is.get();
  } else if (c == 'F') {
      *b = false;
      is.get();
  } else {
    std::cerr << "Read failure in ReadBasicType<bool>, file position is "
              << is.tellg() << ", next char is " << CharToString(c);
  }
}

void ReadBasicType(std::istream &is, bool binary, int *t) {
  assert(t != NULL);
  // Compile time assertion that this is not called with a wrong type.
  if (binary) {
    int len_c_in = is.get();
    if (len_c_in == -1)
      std::cerr << "ReadBasicType: encountered end of stream."<<std::endl;
    char len_c = static_cast<char>(len_c_in), len_c_expected
      = (std::numeric_limits<int>::is_signed ? 1 :  -1)
      * static_cast<char>(sizeof(*t));
    if (len_c !=  len_c_expected) {
      std::cerr << "ReadBasicType: did not get expected integer type, "
                << static_cast<int>(len_c)
                << " vs. " << static_cast<int>(len_c_expected)
                << ".  You can change this code to successfully"
                << " read it later, if needed."<<std::endl;
      // insert code here to read "wrong" type.  Might have a switch statement.
    }
    is.read(reinterpret_cast<char *>(t), sizeof(*t));
  } else {
    if (sizeof(*t) == 1) {
      int16 i;
      is >> i;
      *t = i;
    } else {
      is >> *t;
    }
  }
  if (is.fail()) {
    std::cerr << "Read failure in ReadBasicType, file position is "
              << is.tellg() << ", next char is " << is.peek()<<std::endl;
  }
}
