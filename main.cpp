#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <atomic>
#include "tree.h"

void write_data_to_file(std::unordered_map<std::string,IntDoubleVector>& data,std::string filename) {
    std::ofstream outf{filename};

    if (!outf) {
        std::cerr << "Could not open file "+filename<<'\n';
        return;
    }

    std::vector<std::string> keys{};

    int j=0;
    for (auto const& [key,_]: data) {
        keys.push_back(key);
        if (j==0) outf << key;
        else outf<< ","<<key;
        j++;
    }
    if (keys.size()==0) return;
    outf<<'\n';

    int data_length=get_size(data[keys[0]]);
    for (int i=0; i<data_length; i++) {
        for (int j=0; j<keys.size(); j++) {
            if (std::holds_alternative<std::vector<int>>(data[keys[j]])) {
                int cur_entry{std::get<std::vector<int>>(data[keys[j]])[i]};
                if (j==0) outf<< cur_entry;
                else outf<<","<<cur_entry;
            } else {
                double cur_entry{std::get<std::vector<double>>(data[keys[j]])[i]};
                if (j==0) outf<< cur_entry;
                else outf<<","<<cur_entry;
            }
        }
        if (i!=data_length-1) outf<<'\n';
    }
}



void simulate(Tree& t, double end_time, int seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> uniform_dist(std::nextafter(0.0, 1.0), 1.0);
    double current_time{0};

    std::vector<double> total_rates;
    double total_net_rate{0};
    for (int i=0; i<t.get_number_types(); i++) {
        total_rates.push_back(t.get_population(i) * t.get_rate(i));
        total_net_rate+=t.get_population(i) * t.get_rate(i);
    }

    int type_for_event_to_happen_to{};
    int event_type{};
    double accumulated_event_type_probs{};
    double random_draw{};

    std::vector<std::vector<OffspringPossibility>> offsprings{t.get_offspring_distribution()};

    double next_time_to_print=1;

    TimeSaver time_saver{{1,2,3,4,5,6,7,8,9,10,11,12,13,14}};

    std::unordered_map<std::string,IntDoubleVector> total_population_time{};
    std::unordered_map<std::string,IntDoubleVector> sampled_population_time{};

    while (current_time<=end_time && total_net_rate>0) {
        double dt = -std::log(uniform_dist(rng)) / total_net_rate;
        current_time += dt;
        next_time_to_print-=dt;

        if (next_time_to_print<=0) {
            next_time_to_print=1;
            std::string out_string{'\r'};
            out_string+="Time: "+std::to_string(current_time)+" / "+std::to_string(end_time) + " | ";
            for (int i=0; i<t.get_number_types(); i++) {
                out_string+="T"+std::to_string(i)+": "+std::to_string(t.get_tree_size(i))+" / "+std::to_string(t.get_max_tree(i));
                if (i!=t.get_number_types()-1) out_string+=" | ";
            }
        }
        random_draw = uniform_dist(rng) * total_net_rate;
        total_net_rate=0;
        for (type_for_event_to_happen_to=0; total_net_rate < random_draw; type_for_event_to_happen_to++) {
            total_net_rate+=total_rates[type_for_event_to_happen_to];
        }
        type_for_event_to_happen_to--;

        random_draw=uniform_dist(rng);
        accumulated_event_type_probs=0;
        for (event_type=0; accumulated_event_type_probs < random_draw; event_type++) {
            accumulated_event_type_probs+=offsprings[type_for_event_to_happen_to][event_type].probability;
        }
        event_type--;
        t.replace_cell_random_offspring(type_for_event_to_happen_to,event_type,rng);
        t.perform_save(time_saver,type_for_event_to_happen_to,current_time,event_type);

        while (time_saver.check_last_used_save()) {
            t.get_populations_at_last_save(time_saver,total_population_time);
            t.sample_subtree_binomial({0.5},rng).get_populations_at_last_save(time_saver,sampled_population_time);

            time_saver.increment_last_used_save();
        }

        total_net_rate=0;
        for (int i=0; i<t.get_number_types(); i++) {
            total_rates[i]=t.get_population(i) * t.get_rate(i);
            total_net_rate+=t.get_population(i) * t.get_rate(i);
        }
    }

    write_data_to_file(total_population_time, "total_population/"+std::to_string(seed)+".csv");
    write_data_to_file(sampled_population_time,"sampled_population/"+std::to_string(seed)+".csv");
}

int main() {
    std::vector<OffspringPossibility> type0_offspring{OffspringPossibility{1,{2}}};

    std::vector<double> rates{1};

    constexpr int NUM_SIMULATIONS{1};
    constexpr int NUM_THREADS{1};

    std::atomic<int> completed_simulations{0};
    #pragma omp parallel for num_threads(NUM_THREADS)
    for (int i=0; i<NUM_SIMULATIONS; i++) {
        Tree t{1,0,{type0_offspring},rates,{10'000'000}};
        simulate(t,15,i); 
        int completed=++completed_simulations;
        #pragma omp critical
        {
            int percentage = (100 * completed) / NUM_SIMULATIONS;
            std::cout << "\rProgress: " << percentage << "% (" << completed << "/" << NUM_SIMULATIONS << ") " << std::flush;
        }  
    }
    return 0;
}