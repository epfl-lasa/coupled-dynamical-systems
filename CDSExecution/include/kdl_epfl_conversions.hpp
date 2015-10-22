#ifndef __KDL_EPFL_CONVERSIONS_HPP
#define __KDL_EPFL_CONVERSIONS_HPP

#include <kdl/frames.hpp>
#include <MathLib/Matrix4.h>
#include <MathLib/Matrix3.h>
#include <MathLib/Vector3.h>

KDL::Vector toKDL(const MathLib::Vector3& p);
KDL::Rotation toKDL(const MathLib::Matrix3& M);
KDL::Frame toKDL(const MathLib::Matrix4& T);

MathLib::Vector3 toMathLib(const KDL::Vector& p);
MathLib::Matrix3 toMathLib(const KDL::Rotation& M);
MathLib::Matrix4 toMathLib(const KDL::Frame& T);

#endif  // __KDL_EPFL_CONVERSIONS_HPP
