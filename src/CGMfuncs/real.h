/// @file real.h

#ifndef REAL_H
#define REAL_H

#ifdef _REAL_IS_DOUBLE_
#define real double
#else
#define real float
#endif

extern "C" {
	int sizeof_real_();
}

#endif

