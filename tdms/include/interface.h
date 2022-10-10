/**
 * @file interface.h
 */
#pragma once

class InterfaceComponent{
public:
  bool apply;
  int index;

  InterfaceComponent(const mxArray *ptr, const std::string &name);
};
