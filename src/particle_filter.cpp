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

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
    // Set the number of particles. Initialize all particles to first position (based on estimates of
    // x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    num_particles = 200;

    default_random_engine gen;

    normal_distribution<double> dist_x(x, std[0]);
    normal_distribution<double> dist_y(y, std[1]);
    normal_distribution<double> dist_theta(theta, std[2]);

    Particle temp_particle;

    for(int i = 0; i < num_particles; ++i)
    {
        temp_particle.id     = i;
        temp_particle.x      = dist_x(gen);
        temp_particle.y      = dist_y(gen);
        temp_particle.theta  = dist_theta(gen);
        temp_particle.weight = 1.0;
        temp_particle.associations.clear();
        temp_particle.sense_x.clear();
        temp_particle.sense_y.clear();

        /*
        std::cout << "Particle Id: " << temp_particle.id
                  << ", X: "         << temp_particle.x
                  << ", Y: "         << temp_particle.y
                  << ", Theta: "     << temp_particle.theta
                  << ", Weight: "    << temp_particle.weight
                  << std::endl;
        */

        particles.push_back(temp_particle);

        weights.push_back(temp_particle.weight);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate)
{
    // Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;
    normal_distribution<double> dist_x(0.0,     std_pos[0]);
    normal_distribution<double> dist_y(0.0,     std_pos[1]);
    normal_distribution<double> dist_theta(0.0, std_pos[2]);

    Particle temp_particle;
    double x = 0.0;
    double y = 0.0;
    double theta = 0.0;
    double v_on_omega = 0.0;

    for(unsigned int i = 0; i < particles.size(); ++i)
    {
        temp_particle = particles[i];

        /*
        std::cout << "Original Particle"
                  << ", Id: "    << temp_particle.id
                  << ", X: "     << temp_particle.x
                  << ", Y: "     << temp_particle.y
                  << ", Theta: " << temp_particle.theta
                  << std::endl;
        */

        if(fabs(yaw_rate) > 0.001)
        {
            v_on_omega = velocity / yaw_rate;

            theta = temp_particle.theta + yaw_rate * delta_t;
            x     = temp_particle.x + v_on_omega * (sin(theta) - sin(temp_particle.theta));
            y     = temp_particle.y + v_on_omega * (cos(temp_particle.theta) - cos(theta));

            /*
            std::cout << "Good yaw rate: " << yaw_rate
                      << ", Velocity: "    << velocity
                      << ", Delta_T: "     << delta_t
                      << ", X: "           << x
                      << ", Y: "           << y
                      << ", Theta: "       << theta
                      << std::endl;
            */
        }
        else
        {
            //If dyaw is very small
            theta = temp_particle.theta;
            x     = temp_particle.x + velocity * delta_t * cos(theta);
            y     = temp_particle.y + velocity * delta_t * sin(theta);

            /*
            std::cout << "Small yaw rate: " << yaw_rate
                      << ", Velocity: "     << velocity
                      << ", Delta_T: "      << delta_t
                      << ", X: "            << x
                      << ", Y: "            << y
                      << ", Theta: "        << theta
                      << std::endl;
            */
        }

        //Adding Gaussian noise for X, Y and Theta
        particles[i].x     = x + dist_x(gen);
        particles[i].y     = y + dist_y(gen);
        particles[i].theta = theta + dist_theta(gen);

        /*
        std::cout << "Predicted Particle"
                  << ", Id: "    << particles[i].id
                  << ", X: "     << particles[i].x
                  << ", Y: "     << particles[i].y
                  << ", Theta: " << particles[i].theta
                  << std::endl;
        */
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations)
{
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks)
{
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

    /*
    static unsigned int counter = 0;

        for(unsigned int j = 0; j < observations.size(); ++j)
        {
            std::cout << "Counter: " << (++counter) << ", Observation " << observations[j].id << ", X: "  << observations[j].x << ", Y: "  << observations[j].y << std::endl;
        }
    */

    const double GAUSSIAN_NORM = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
    const double ONE_ON_2_SIGMAX_2   = 1.0 / (2.0 * std_landmark[0] * std_landmark[0]);
    const double ONE_ON_2_SIGMAY_2   = 1.0 / (2.0 * std_landmark[1] * std_landmark[1]);

    for(unsigned int i = 0; i < particles.size(); ++i)
    {
        std::cout << "Particle " << i << ", X: " << particles[i].x << ", Y: " << particles[i].y << ", Theta: " << particles[i].theta << std::endl;

        double final_weight = 1.0;

        for(unsigned int j = 0; j < observations.size(); ++j)
        {
            std::cout << "Observation " << observations[j].id << ", X: "  << observations[j].x << ", Y: "  << observations[j].y << std::endl;

            double x_obs_map = particles[i].x + cos(particles[i].theta) * observations[j].x - sin(particles[i].theta) * observations[j].y;
            double y_obs_map = particles[i].y + sin(particles[i].theta) * observations[j].x + cos(particles[i].theta) * observations[j].y;

            std::cout << "Converted observation, X: " << x_obs_map << ", Y: " << y_obs_map << std::endl;

            //Reset nearest distance to uninitialised and nearest landmark not found
            double temp_dist = 0.0;
            double nearest_dist = -1.0;
            //double found_nearest_landmark = false;

            double x_nearest_landmark = 0.0;
            double y_nearest_landmark = 0.0;

            for(unsigned int k = 0; k < map_landmarks.landmark_list.size(); ++k)
            {
                temp_dist = dist(x_obs_map, y_obs_map, map_landmarks.landmark_list[k].x_f, map_landmarks.landmark_list[k].y_f);

                //if( !(temp_dist > sensor_range) )
                //{
                    if( (nearest_dist == -1.0) || (temp_dist < nearest_dist) )
                    {
                        nearest_dist = temp_dist;
                        x_nearest_landmark = map_landmarks.landmark_list[k].x_f;
                        y_nearest_landmark = map_landmarks.landmark_list[k].y_f;
                        //found_nearest_landmark = true;
                    }
                //}
                //else
                //{
                //    std::cout << "Landmark "   << k << " (" << map_landmarks.landmark_list[k].x_f << ", " << map_landmarks.landmark_list[k].y_f
                //               << "), has distance " << temp_dist << " > sensor range " << sensor_range << "! Ignored!" << std::endl;
                //}
            }

            //if(found_nearest_landmark)
            //{
                double delta_x    = (x_obs_map - x_nearest_landmark);
                double delta_y    = (y_obs_map - y_nearest_landmark);
                double obs_weight = GAUSSIAN_NORM * exp( -(ONE_ON_2_SIGMAX_2 * delta_x * delta_x + ONE_ON_2_SIGMAY_2 * delta_y * delta_y));

                //if(0.0 == final_weight)
                //{
                //    final_weight = obs_weight;
                //}
                //else
                //{
                    final_weight *= obs_weight;
                //}

                std::cout << "Weight from observation: " << obs_weight << ", Updated final weight: " << final_weight << std::endl;
            //}
            //else
            //{
                //std::cout << "Nearest landmark within sensor range " << sensor_range << " not found!" << std::endl;
            //}
        }

        particles[i].weight = final_weight;

        std::cout << "Particle " << i << " final weight: " << particles[i].weight << std::endl;
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

void ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
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
