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
    num_particles = 15;

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
    default_random_engine gen;
    normal_distribution<double> dist_x(0.0, std_pos[0]);
    normal_distribution<double> dist_y(0.0, std_pos[1]);
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
        particles[i].x     = x     + dist_x(gen);
        particles[i].y     = y     + dist_y(gen);
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
    //Not implemented
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks)
{
    const double GAUSSIAN_NORM       = 0.5 / (M_PI * std_landmark[0] * std_landmark[1]);
    const double ONE_ON_2_SIGMAX_2   = 0.5 / (std_landmark[0] * std_landmark[0]);
    const double ONE_ON_2_SIGMAY_2   = 0.5 / (std_landmark[1] * std_landmark[1]);

    //Find nearest landmarks within sensor range around each particle
    std::vector< std::vector<Map::single_landmark_s> > nearest_landmarks;

    for(unsigned int pt_id = 0; pt_id < particles.size(); ++pt_id)
    {
        std::vector<Map::single_landmark_s>landmark_vector;

        for(unsigned int lm_id = 0; lm_id < map_landmarks.landmark_list.size(); ++lm_id)
        {
            double particle_dist = dist(particles[pt_id].x, particles[pt_id].y, (map_landmarks.landmark_list[lm_id]).x_f, (map_landmarks.landmark_list[lm_id]).y_f);

            if(particle_dist < sensor_range)
            {
                landmark_vector.push_back(map_landmarks.landmark_list[lm_id]);
            }
        }

        nearest_landmarks.push_back(landmark_vector);
    }

    //Vector of multivariate Gaussian probabilities for each particle
    std::vector<double> gaussian_probs;
    double sum_particle_raw_weight = 0.0;

    for(unsigned int pt_id = 0; pt_id < particles.size(); ++pt_id)
    {
        /*
        std::cout << "Update weight for Particle[" << pt_id
                  << "]: (" << particles[pt_id].x
                  << ", "   << particles[pt_id].y
                  << ", "   << particles[pt_id].theta
                  << ", "   << particles[pt_id].weight
                  << ")"    << std::endl;
        */

        //Reset vector of Gaussian probabilities
        gaussian_probs.clear();

        for(unsigned int obs_id = 0; obs_id < observations.size(); ++obs_id)
        {
            //std::cout << " - Observation[" << obs_id << "]: ("  << observations[obs_id].x << ", "  << observations[obs_id].y << ")" << std::endl;

            //Convert observation to map coordinates w.r.t the current particle
            double x_obs_map = particles[pt_id].x + cos(particles[pt_id].theta) * observations[obs_id].x - sin(particles[pt_id].theta) * observations[obs_id].y;
            double y_obs_map = particles[pt_id].y + sin(particles[pt_id].theta) * observations[obs_id].x + cos(particles[pt_id].theta) * observations[obs_id].y;

            //std::cout << " - Converted to Map coordinates (" << x_obs_map << ", " << y_obs_map << ")" << std::endl;

            //Search for nearest landmark around the current particle
            double dist_obs              =   0.0;
            double nearest_dist          =  -1.0;
            double x_nearest_landmark    =   0.0;
            double y_nearest_landmark    =   0.0;
            bool found_nearest_landmark  = false;

            for(unsigned int lm_id = 0; lm_id < (nearest_landmarks[pt_id].size()); ++lm_id)
            {
                dist_obs = dist(x_obs_map, y_obs_map, (nearest_landmarks[pt_id][lm_id]).x_f, (nearest_landmarks[pt_id][lm_id]).y_f);

                if( (-1.0 == nearest_dist) || (dist_obs < nearest_dist) )
                {
                    nearest_dist           = dist_obs;
                    x_nearest_landmark     = (nearest_landmarks[pt_id][lm_id]).x_f;
                    y_nearest_landmark     = (nearest_landmarks[pt_id][lm_id]).y_f;
                    found_nearest_landmark = true;
                }
            }

            //Compute a multivariate Gaussian probability for the current observation w.r.t the current particle
            if(found_nearest_landmark)
            {
                double delta_x    = (x_obs_map - x_nearest_landmark);
                double delta_y    = (y_obs_map - y_nearest_landmark);
                double obs_weight = GAUSSIAN_NORM * exp( -((ONE_ON_2_SIGMAX_2 * delta_x * delta_x) + (ONE_ON_2_SIGMAY_2 * delta_y * delta_y)) );

                gaussian_probs.push_back(obs_weight);

                /*
                std::cout << " - Nearest Landmark"
                          << " (" << x_nearest_landmark << ", " << y_nearest_landmark << ")"
                          << ", Raw weight from Observation: "  << obs_weight << std::endl;
                */
            }
        }

        //Multiply multivariate Gaussian probabilities for observations w.r.t the current particle.
        if(gaussian_probs.size())
        {
            particles[pt_id].weight = 1.0;

            for(unsigned int gp_id = 0; gp_id < gaussian_probs.size(); ++gp_id)
            {
                particles[pt_id].weight *= gaussian_probs[gp_id];
            }
        }
        else
        {
            particles[pt_id].weight = 0.0;
        }

        //Sum of all raw particle weights
        sum_particle_raw_weight += particles[pt_id].weight;
    }

    //Normalise particle weights
    if(sum_particle_raw_weight > 0.0)
    {
        for(unsigned int pt_id = 0; pt_id < particles.size(); ++pt_id)
        {
            particles[pt_id].weight /= sum_particle_raw_weight;
            weights[pt_id]           = particles[pt_id].weight;

            //std::cout << " - Normalised weight of Particle[" << pt_id<< "]: " << particles[pt_id].weight << std::endl;
        }
    }
}

void ParticleFilter::resample()
{
    //C++ re-implementation of Python Resampling wheel from the lesson.

    default_random_engine gen;

    uniform_int_distribution<int> uniform_random_id(0, num_particles - 1);

    double max_weight = *(std::max_element(weights.begin(), weights.end()));
    uniform_real_distribution<double> uniform_random_weight(0.0, (2.0 * max_weight));

    int index = uniform_random_id(gen);
    double beta = 0.0;

    vector<Particle> temp_particles;

    for(unsigned int pt_id = 0; pt_id < particles.size(); ++pt_id)
    {
        beta += uniform_random_weight(gen);

        while(beta > weights[index])
        {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }

        temp_particles.push_back(particles[index]);
    }

    particles = temp_particles;
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

