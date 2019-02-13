/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <map>

#include "helper_functions.h"


using std::string;
using std::vector;
using namespace std;

default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  num_particles = 100;  // TODO: Set the number of particles
  //creates a normal (Gaussian) distribution for x
  normal_distribution<double> dist_x(x,std[0]);
  //creates a normal (Gaussian) distribution for y 
  normal_distribution<double> dist_y(y,std[1]);
  //creates a normal (Gaussian) distribution for theta
  normal_distribution<double> dist_theta(theta,std[2]);
  
  for (int i = 0; i < num_particles; ++i) {
    Particle particle;
    particle.id = i;
    particle.x = dist_x(gen);
    particle.y = dist_y(gen);
    particle.theta = dist_theta(gen);
    particle.weight = 1.0;
    particles.push_back(particle);
    weights.push_back(particles[i].weight);
  }
  is_initialized = true;
}
void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  normal_distribution<double> noise_x(0,std_pos[0]);
    //creates a normal (Gaussian) distribution for y 
  normal_distribution<double> noise_y(0,std_pos[1]);
    //creates a normal (Gaussian) distribution for theta
  normal_distribution<double> noise_theta(0,std_pos[2]);
  
  double xf, yf, thetaf;

  for (int i = 0; i < num_particles; i++) {
    double thetaf = particles[i].theta;
    if (fabs(yaw_rate)< 0.00001) {
      xf = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      yf = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      
    } else {
      xf = particles[i].x + (velocity/yaw_rate)*(sin(particles[i].theta+delta_t*yaw_rate)-sin(particles[i].theta));
      yf = particles[i].y + (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+delta_t*yaw_rate));
      thetaf = particles[i].theta + yaw_rate*delta_t;
    }
      //creates a normal (Gaussian) distribution for x
      particles[i].x = xf + noise_x(gen);
      particles[i].y = yf + noise_y(gen);
      particles[i].theta = thetaf + noise_theta(gen);
      // std::cout << "x:" << particles[i].x << "\n" << "y:" << particles[i].y << "\n";
  } 
}

  

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (unsigned int i = 0; i < observations.size(); i++) {
    double min_dist = 9999999.;
    LandmarkObs obs = observations[i];
    double o_x = obs.x;
    double o_y = obs.y;
    for (unsigned int j = 0;j < predicted.size();j++) {
      LandmarkObs prd = predicted[j];
      double p_x = prd.x;
      double p_y = prd.y;
      if (dist(o_x,o_y,p_x,p_y) < min_dist) {
	min_dist = dist(o_x,o_y,p_x,p_y);
	observations[i].id = prd.id;
      }
    }   
  }
}
// create a function to calculate the multi_gauss probability


double multiv_prob(double sig_x, double sig_y,double x_obs, double y_obs,double mu_x,double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);
  
  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x,2) / (2 * pow(sig_x, 2))) + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
  
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
  
  return weight;
}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  double weights_total = 0.0; 
  vector<Map::single_landmark_s> ldmks = map_landmarks.landmark_list;
  for (unsigned int i = 0; i < particles.size(); i++) {
    vector<LandmarkObs> predicted;
    double pr_x = particles[i].x;
    double pr_y = particles[i].y;
    double pr_theta = particles[i].theta;
    //create a vector of LandmarkObs contains the predicted landmarks within the lidar sensor range
    for (unsigned int j = 0; j < ldmks.size(); j++) {
      LandmarkObs pre_ldm;
      double ld_x = ldmks[j].x_f;
      double ld_y = ldmks[j].y_f;
      if (dist(pr_x,pr_y,ld_x,ld_y) < sensor_range) {
	pre_ldm.id = ldmks[j].id_i;
	pre_ldm.x = ld_x;
	pre_ldm.y = ld_y;
	predicted.push_back(pre_ldm);
      } 
    }
    //create a vector of Landmarks contains the observations landmarks transformed in the map coordinates
    vector<LandmarkObs> tf_observations;
    for (unsigned int k = 0; k < observations.size(); k++) {
      LandmarkObs tf_ob;
      tf_ob.x = pr_x + cos(pr_theta)*observations[k].x - sin(pr_theta)*observations[k].y;
      tf_ob.y = pr_y + sin(pr_theta)*observations[k].x + cos(pr_theta)*observations[k].y;
      tf_ob.id = observations[k].id;
      tf_observations.push_back(tf_ob);
    }
    // association the predicted landmarks with the measurements
    dataAssociation(predicted, tf_observations);
    
    // update the weight of each particle based on the mutivariate Gaussian probability function
    for (unsigned int m = 0; m < tf_observations.size(); m++) {
      for (unsigned int n = 0; n < predicted.size(); n++) {
	if (tf_observations[m].id == predicted[n].id) {
	  double weight = multiv_prob(std_landmark[0],std_landmark[1],tf_observations[m].x,tf_observations[m].y,predicted[n].x,predicted[n].y);
	  if (weight == 0) {
	    particles[i].weight *= 0.00001;
	  } else {
	    particles[i].weight *= weight;
	  }
	  //std::cout << "weight:" << weight;
	}
      }
    }
    weights_total += particles[i].weight;
    
  }
  //normalize the weights vector
  for (unsigned int l = 0; l < particles.size(); l++) {
    particles[l].weight /= weights_total;
    weights[l] =  particles[l].weight;
  }
  std::cout << "weights_total:" << weights_total << "\n";
}
//cout << weights;
void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> rsp_particles;
  std::default_random_engine gen;
  std::discrete_distribution<int> d(weights.begin(), weights.end());
  for(unsigned int i=0; i<particles.size(); i++) {
    int index = d(gen);
    rsp_particles.push_back(particles[index]);
  }
  particles = rsp_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}