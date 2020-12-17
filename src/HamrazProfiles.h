#ifndef HAMRAZPROFILE_H
#define HAMRAZPROFILE_H

#include <Rcpp.h>
#include <lidR/Point.h>

namespace Hamraz
{
  class Profile
  {
    public:
      Profile(std::vector<lidR::PointXYZR*>& Data, lidR::PointXYZR Center, double Angle, double Radius, double Width, int Sensitivity, double MDCW, double Epsilon, double CLc, double CLs, double Oc, double Os);
      ~Profile();
      Rcpp::List to_R();

    public:
      double angle;
      double width;
      int sensitivity;
      double mdcw;
      double epsilon;
      double clc;
      double cls;
      double oc;
      double os;
      lidR::PointXYZR center;
      std::vector<lidR::PointXYZR*> points;
      std::vector<lidR::PointXYZR*> points_no_gaps;
      std::vector<lidR::PointXYZR*> points_no_boundaries;
      lidR::PointXYZR extremityPoint;
      std::vector<int> localMinimaIndex;

    private:
      void extract_profile(std::vector<lidR::PointXYZR*>& Data);
      void find_inter_tree_gaps();
      void find_boundary();
      void find_local_minima();
      double IQR(std::vector<double>);
      double median(std::vector<double>);
      double steepness(std::vector<lidR::PointXYZR*> &subProfile);
      void extract_points_prior(std::vector<lidR::PointXYZR*> &subProfile, double limit, std::vector<lidR::PointXYZR*> &subProfileSubset);

  };

  bool operator<(Profile const& a, Profile const& b);

  class ProfilesManager
  {
    public:
      ProfilesManager(std::vector<lidR::PointXYZR*>& data, lidR::PointXYZR Center, double Radius, double Width, int Sensitivity, double MDCW, double Epsilon, double CLc, double CLs, double Oc, double Os);
      ~ProfilesManager();
      void add_next_profiles(std::vector<lidR::PointXYZR*>& data);
      std::vector<lidR::PointXYZ> get_polygon();
      Rcpp::List to_R();

    public:
      double chord;

    private:
      int sensitivity;
      double alpha;
      double rmax;
      double radius;
      double width;
      double mdcw;
      double epsilon;
      double clc;
      double cls;
      double oc;
      double os;
      lidR::PointXYZR center;
      std::vector<Profile> profiles;
      std::vector<int> localMinimaIndex;
  };
}

#endif //HAMRAZPROFILE_H
