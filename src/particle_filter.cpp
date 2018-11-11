/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// This line creates a normal (Gaussian) distribution for x
	num_particles = 100;
	is_initialized = true;
	// default_random_engine gen;
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	Particle p;
	for (int i = 0; i < num_particles; i++) {
		 p.id = i;
		 p.x = dist_x(gen);
		 p.y = dist_y(gen);
		 p.theta = dist_theta(gen);	 
		 p.weight = 1;
		 weights.push_back(p.weight);
		 particles.push_back(p);
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	// default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);
	// store variable to use for all particles
	double vel_by_yaw_rate = velocity/yaw_rate;

	if(fabs(yaw_rate)>0.0001) {
		for(int i = 0; i < num_particles; i++) {
			particles[i].x += vel_by_yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) -  sin(particles[i].theta)) + dist_x(gen);
			particles[i].y += vel_by_yaw_rate * (cos(particles[i].theta) -  cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
			particles[i].theta += yaw_rate*delta_t + dist_theta(gen);
		}
	}
	else {
		for(int i = 0; i < num_particles; i++) {
			particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
      		particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
      		particles[i].theta += dist_theta(gen);
		}
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for(int i = 0; i < observations.size(); i++) {
		double min_dist = numeric_limits<double>::max();
		double curr_dist; // to store distance values in each iteration
		for(int j = 0; j < predicted.size(); j++) {
			curr_dist = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
			if(curr_dist < min_dist) {
				/* note that the id here is just the index in predicted, otherwise to find the closest
				landmark, will have to loop through predicted vector again looking for the landmark ID */
				observations[i].id = j;
				min_dist = curr_dist;
			}
		}
	}


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	// loop over all particles
	// cout << "updating weights" << "\n";
	for(int i = 0; i < num_particles; i++) {
		// variables for storing landmark id, x and y
		int landmark_id;
		float landmark_x, landmark_y;
		// vector to store landmarks within distance sensor_range of particle
		vector<LandmarkObs> predicted;
		// looping over all the landmarks to collect landmarks within sensor range of particle
		for(int j = 0; j < map_landmarks.landmark_list.size(); j++) {
			landmark_id = map_landmarks.landmark_list[j].id_i;
			landmark_x = map_landmarks.landmark_list[j].x_f;
			landmark_y = map_landmarks.landmark_list[j].y_f;
			if(dist(landmark_x, landmark_y, particles[i].x, particles[i].y)<= sensor_range) 
				predicted.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
		}

		// transformed coordinates vector
	    vector<LandmarkObs> transformed_observations;
	    // coordinates in map system
	    double x_m, y_m;
	    // transform all observations to the map coordinate system
	    for(int j = 0; j < observations.size(); j++) {
	      x_m = cos(particles[i].theta)*observations[j].x - sin(particles[i].theta)*observations[j].y + particles[i].x;
	      y_m = sin(particles[i].theta)*observations[j].x + cos(particles[i].theta)*observations[j].y + particles[i].y;
	      transformed_observations.push_back(LandmarkObs{observations[j].id, x_m, y_m});
	    }

	    // perform data association
	    dataAssociation(predicted, transformed_observations);

	    double weight = 1;
	    int index_landmark;
	    // finally update the weights using the Gaussian assumption
	    // if no landmark is within sensor range, that particle is useless
	    if(predicted.size() == 0) {
	    	weight = 0;
	    }
	    else {
	    	for(int j = 0; j < transformed_observations.size(); j++) {
		    	index_landmark = transformed_observations[j].id;
		    	landmark_x = predicted[index_landmark].x;
		    	landmark_y = predicted[index_landmark].y;
		    	weight *= eval_gaussian(landmark_x, landmark_y, std_landmark[0], std_landmark[1], 
		    								transformed_observations[j].x, transformed_observations[j].y);
	    	}
	    
	    }
	
	    particles[i].weight = weight;
	    weights[i] = weight;

	}
}
	



void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	// to store new particles
	// default_random_engine gen;

	// cout << "resampling" << "\n";

	vector<Particle> new_particles;
	double max_weight = numeric_limits<double>::min();
	// loop for finding max weight
	for (int i = 0; i < num_particles; i++) {
		max_weight = max(max_weight, particles[i].weight);
  	}
	double beta = 0;
	// for uniform integer values
	uniform_int_distribution<int> unif_int(0, num_particles-1);

	// for uniform real distributions
	uniform_real_distribution<double> unif_real(0.0,1.0);

	// generate a random index
	int index = unif_int(gen);
	// re-sampling wheel implementation
	for(int i = 0; i < num_particles; i++) {
		beta = beta + unif_real(gen) * 2 * max_weight;
		while(particles[index].weight < beta) {
			beta = beta - particles[index].weight;
			index = (index + 1) % num_particles;
		}
		new_particles.push_back(particles[index]);
	}
	particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
