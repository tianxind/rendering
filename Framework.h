#ifndef FRAMEWORK_H
#define FRAMEWORK_H

// GLEW must be included first, if we need it.
#ifdef _WIN32
//#define GLEW_STATIC
#define FRAMEWORK_USE_GLEW
#include "glew.h"
#endif

#ifdef __linux__
#define FRAMEWORK_USE_GLEW
#include "glew.h"
#endif

#include <memory>
#include <iostream>

#endif
