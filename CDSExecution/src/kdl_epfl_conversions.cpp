#include <kdl_epfl_conversions.hpp>

KDL::Vector toKDL(const MathLib::Vector3& p)
{
  return KDL::Vector(p.cx(), p.cy(), p.cz());
}

KDL::Rotation toKDL(const MathLib::Matrix3& M)
{
  return KDL::Rotation(M.At(0,0), M.At(0,1), M.At(0,2), M.At(1,0), M.At(1,1), M.At(1,2), M.At(2,0), M.At(2,1), M.At(2,2));
}

KDL::Frame toKDL(const MathLib::Matrix4& T)
{
  return KDL::Frame(toKDL(T.GetOrientation()), toKDL(T.GetTranslation()));
}

MathLib::Vector3 toMathLib(const KDL::Vector& p)
{
  return MathLib::Vector3(p.x(), p.y(), p.z());
}

MathLib::Matrix3 toMathLib(const KDL::Rotation& M)
{
  return MathLib::Matrix3(M.data);
}

MathLib::Matrix4 toMathLib(const KDL::Frame& T)
{
  return MathLib::Matrix4(toMathLib(T.M), toMathLib(T.p));
}
