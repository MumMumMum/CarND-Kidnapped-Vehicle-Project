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
//static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).

    if(!is_initialized ) {
        num_particles = 9;

        const double std_x = std[0];
        const double std_y = std[1];
        double std_theta = std[2];
        particles.clear();
        weights.clear();
        default_random_engine gen;
        normal_distribution<double> dist_x(x, std_x);
        normal_distribution<double> dist_y(y, std_y);
        normal_distribution<double> dist_theta(theta, std_theta);
        //normal_distribution<double > dist_noise(0,0.01);
        //double newx,newy,newtheta;
        for (int i = 0; i < num_particles; i++) {
            Particle p;

            p.x = dist_x(gen);
            p.y = dist_y(gen);
            p.theta = dist_theta(gen);
            p.weight = 1;
            particles.push_back(p);
            //printf("the initalized particles are x y theta, %f %f %f ",p.x,p.y,p.theta);
            //cout<<endl;
            weights.push_back(1.0);
        }
    }
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    default_random_engine gen;
    const double std_x = std_pos[0];
    const double std_y = std_pos[1];
    const double std_theta = std_pos[2];
    double x,y,theta;
    if (fabs(yaw_rate) < 0.00001) {
        for (int i = 0; i < num_particles; i++) {
            x = particles[i].x;
            y = particles[i].y;
            theta = particles[i].theta;
            x = x + velocity * delta_t * cos(theta);
            y = y + velocity * delta_t * sin(theta);
            theta = theta + (yaw_rate * delta_t);

            normal_distribution<double> dist_x(0, std_x);
            normal_distribution<double> dist_y(0, std_y);
            normal_distribution<double> dist_theta(0, std_theta);
            //add random noise
            particles[i].x     = dist_x(gen) + x;
            particles[i].y     = dist_y(gen) + y;
            particles[i].theta = dist_theta(gen)+theta;

        }
    }
    else{

        for (int i = 0; i < num_particles; i++) {
            x = particles[i].x;
            y = particles[i].y;
            theta = particles[i].theta;
            x = x +
                (velocity / yaw_rate) * (sin(theta + (yaw_rate * delta_t)) - sin(theta));
            y = y +
                (velocity / yaw_rate) * (cos(theta) - cos(theta + (yaw_rate * delta_t)));
            theta = theta + (yaw_rate * delta_t);

            normal_distribution<double> dist_x(0, std_x);
            normal_distribution<double> dist_y(0, std_y);
            normal_distribution<double> dist_theta(0 , std_theta);
            //add random noise
            particles[i].x     = dist_x(gen) + x;
            particles[i].y     = dist_y(gen) + y;
            particles[i].theta = dist_theta(gen)+theta;
            }

        }

        //std::cout<<"---------------------Prediction Calc------------------------"<<std::endl;
        //std::cout << "Velocity:"<<velocity<<";Yawrate:"<<yaw_rate<<std::endl;
        //std ::cout <<"Predicted:x,y,theta:"<<particles[i].x<<" , "<<particles[i].y<<" , "<<particles[i].theta<<std::endl;
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs> &observations) {

    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.
    for (int i = 0; i < observations.size(); i++) {
        int map_id = 0;
        LandmarkObs obsv = observations[i];
        //double min_dist = std::numeric_limits<double>::max();
        double min_dist = 1000.00;
        for (int j = 0; j < predicted.size(); j++) {

            LandmarkObs pred = predicted[j];
            double meas_dist = dist(obsv.x, obsv.y, pred.x, pred.y);
            cout<<"meas_dist:  "<<meas_dist<<" , "<< "LandMark(x,y): "<<pred.x<<","<<pred.y<<", "<<pred.id<<endl;
            if (meas_dist < min_dist) {
                min_dist = meas_dist;
                //cout<<"meas_dist:  "<<meas_dist<<" , "<< "LandMark(x,y): "<<pred.x<<","<<pred.y<<", "<<pred.id<<endl;
                map_id = pred.id;
            }


        }

        observations[i].id = map_id;
        std::cout <<"Landmark Index: "<<observations[i].id <<
        " TObservation(x,y): "<<"("<<observations[i].x<<","<<observations[i].y<<")"<<std::endl;

    }


}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   std::vector<LandmarkObs> observations, Map map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d. <  <<
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html

    std::vector<LandmarkObs> landMarks_inRange;
    std::vector<LandmarkObs> trans_maped;
    double total_weight;
    const double s_x = std_landmark[0];
    const double s_y = std_landmark[1];
    const double c1 = (1 / (2 * M_PI * s_x * s_y));
    const double c2 = 2 * (s_x * s_x);//(2 * pow(s_x, 2));
    const double c3 = 2 * (s_y * s_y);//(2 * pow(s_y, 2));
    for (int i = 0; i < num_particles; i++) {

        Particle p = particles[i];
        double p_x = p.x;
        double p_y = p.y;
        double p_theta = p.theta;
        // so find all the landmark which are in range of sensor
        std::cout << "--------------------LandMark in range------------------------" << std::endl;
        landMarks_inRange.clear();
        for (int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            double lx = map_landmarks.landmark_list[j].x_f;
            double ly = map_landmarks.landmark_list[j].y_f;
            int    lid = map_landmarks.landmark_list[j].id_i;
            double landMark_dist = dist(lx, ly, p.x, p.y);
            if (landMark_dist < sensor_range) {

                landMarks_inRange.push_back(LandmarkObs{lid, lx, ly});
                cout << "Land Mark in Range:id,x,y:" << lid << "," << lx << "," << ly << endl;
           }

        }
        std::cout << "--------------------Transformed Obsv------------------------" << std::endl;
        //mark the particle in the cordinate of the map
        trans_maped.clear();
        for (int j = 0; j < observations.size(); j++) {
            double x = p_x + observations[j].x * cos(p_theta) - observations[j].y * sin(p_theta);
            double y = p_y + observations[j].x * sin(p_theta) + observations[j].y * cos(p_theta);
            trans_maped.push_back(LandmarkObs{observations[j].id, x, y});
            std::cout << "Obs(x,y)" << "(" << observations[j].x << "," << observations[j].y << ")";
            std::cout << "-->TObs(x,y)" << "(" << x << "," << trans_maped[j].y << ")" << std::endl;
        }
        std::cout << "---------------------Assosciations------------------------" << std::endl;
        dataAssociation(landMarks_inRange, trans_maped);


        std::cout << "---------------------Weights Calc------------------------" << std::endl;
        double total_prob = 1.0;
        for (unsigned int j = 0; j < trans_maped.size(); j++) {

            // placeholders for observation and associated prediction coordinates
            // total_prob= 1.0;
            double o_x    = trans_maped[j].x;
            double o_y    = trans_maped[j].y;
            int    o_id   = trans_maped[j].id;
            double pr_x, pr_y;

            // get the x,y coordinates of the prediction associated with the current observation
            for (unsigned int k = 0; k < landMarks_inRange.size(); k++) {
                if (landMarks_inRange[k].id == o_id) {
                    pr_x = landMarks_inRange[k].x;
                    pr_y = landMarks_inRange[k].y;
                    //cout<<"id in "<< trans_maped[j].id<<endl;
                    std::cout << "LandmarkIndex :" << landMarks_inRange[k].id << std::endl;
                }
            }

            std::cout << "Landmark(x,y):" << "(" << pr_x << "," << pr_y << "); Particle(x,y):(" << o_x << "," << o_y
                      << ")" << std::endl;

            //std::cout << "the pr_arr and pr x ny is"<<pr_x,pr_y,pr_arr[0],pr_arr[1]<<endl

            // calculate weight for this observation with multivariate Gaussian
            double dx = o_x - pr_x;
            double dy = o_y - pr_y;
            double dxdx = (dx * dx) / c2;
            double dydy = (dy * dy) / c3;
            double c4 = exp(-(dxdx + dydy));
            double obs_w = c1 * c4;
            //cout<<"total_prob befor:"<<total_prob<<"for observation :"<<j<<endl;
            total_prob *= obs_w;
            std::cout << "dx:" << dx << ";" << "dy:" << dy << " prob for observation:"<<obs_w<<"total Prob: "<<total_prob<<std::endl;

        }//j ends

        particles[i].weight = total_prob;
        weights.push_back(total_prob);
        std::cout << i << " :particles[i].weight " << "=" << particles[i].weight << "calc total prob = " << total_prob
                  << std::endl;
    }//i ends
    //normalize the wt
    //total_weight = 1.0;
    total_weight = std::accumulate(weights.begin(), weights.end(), 0.0);
    //weights.clear();
    for (int i = 0; i < num_particles; i++) {
        particles[i].weight = particles[i].weight / total_weight;
        weights[i] = particles[i].weight;
        //std::cout <<"Weight " <<i<< "="<< weights[i]<<std::endl;

    }
}//fun ends

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    double max_weight = *max_element(weights.begin(), weights.end());

    //std::uniform_real_distribution<double> distr_weight(0.0, max_weight);
    default_random_engine gen;
    std::discrete_distribution< > d(weights.begin(), weights.end());
    vector<Particle> particle_sampled;
    for (int i = 0; i < num_particles; i++) {
        auto index = d(gen);
        index = (index) % num_particles;
        particle_sampled.push_back(particles[index]);
        //cout<<"get wt done"<<index<<endl;
    }
    particles.clear();
    particles = std::move(particle_sampled);
    std::cout << "---------------------resmaple Calc  done------------------------" << std::endl;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x,
                                         std::vector<double> sense_y) {
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();

    particle.associations = associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best) {
    vector<int> v = best.associations;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseX(Particle best) {
    vector<double> v = best.sense_x;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}

string ParticleFilter::getSenseY(Particle best) {
    vector<double> v = best.sense_y;
    stringstream ss;
    copy(v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length() - 1);  // get rid of the trailing space
    return s;
}