/*
 * engine.h
 *
 *  Created on: Sep 19, 2017
 *      Author: tao
 */

#ifndef ENGINE_H_
#define ENGINE_H_

#include "engine/vpcontainer.h"
#include "util/data-interface.h"
#include "base/common.h"
#include "feature/mfcc.h"
#include "feature/cmvn.h"
#include <map>

class VPengine{

public:

	VPengine():container_(NULL),datainterface_(NULL),mfcc_(NULL),cmn_(NULL),threshold_(-4.46516){};

	~VPengine();

	void release();

	void initialize(VPcontainer* container);

	void sign(std::string key, int32 length);

	void check(int32 length);

	void put_data(short* short_data , int32 length);

	void data_end();

private:

	VPcontainer* container_ ;
	DataInterface* datainterface_;
	Mfcc* mfcc_ ;
	CMN* cmn_ ;
	std::map<std::string, Vector> database_ ;
	std::map<std::string, int32> count_ ;
	BaseFloat threshold_;

};



#endif /* ENGINE_H_ */
