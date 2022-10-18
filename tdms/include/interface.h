/**
 * @file interface.h
 */
#pragma once

#include <string>

#include "mat_io.h"

class InterfaceComponent{
public:
  bool apply;
  int index;

  InterfaceComponent(const mxArray *ptr, const std::string &name);
};
