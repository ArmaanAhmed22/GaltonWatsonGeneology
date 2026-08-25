#ifndef TREE_H
#define TREE_H
#include <random>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <variant>
#include <map>
#include <string>
#include <cmath>

using IntDoubleVector = std::variant<std::vector<int>,std::vector<double>>;

int get_size(IntDoubleVector& arr) {
    if (std::holds_alternative<std::vector<int>>(arr)) {
        return std::get<std::vector<int>>(arr).size();
    } else {
        return std::get<std::vector<double>>(arr).size();
    }
}



struct OffspringPossibility {
    double probability;
    std::vector<int> type_to_number;
};

struct SaveSnapshot {
    const std::vector<int> alive_pointer_left;
    const std::vector<int> alive_pointer_right;
    const double time;
    const short type_of_cell_of_interest;
    const int id_of_cell_of_interest;

    int get_total_alive() {
        int tots{};
        for (int i=0; i<alive_pointer_left.size(); i++) {
            tots+=alive_pointer_right[i]-alive_pointer_left[i];
        }
        return tots;
    }
};

struct SaveDataSnapshot {
    std::vector<std::vector<int>> id_to_id_ancestry{};
    std::vector<std::vector<short>> id_to_ancestor_type{};
};

template <typename T>
class Saver;

class Tree {
    private:
        const std::vector<int> MAX_TREE_SIZE_BY_TYPE;
        const short STARTING_TYPE;
        std::vector<std::vector<int>> id_to_id_ancestry{};
        std::vector<std::vector<int>> id_to_founder{};
        std::vector<std::vector<short>> id_to_ancestor_type{};
        std::vector<std::vector<double>> id_to_birth_time{};
        std::vector<std::vector<double>> id_to_death_time{};
        int number_types;
        const std::vector<std::vector<OffspringPossibility>> offspring_distribution;
        std::vector<double> rates;

        std::vector<int> alive_pointer_left{};
        std::vector<int> alive_pointer_right{};
    public:
        Tree(int number_types, short starting_type, std::vector<std::vector<OffspringPossibility>> offspring_distribution, std::vector<double> rates, std::vector<int> max_tree_size_by_type): number_types{number_types},MAX_TREE_SIZE_BY_TYPE{max_tree_size_by_type}, offspring_distribution{offspring_distribution}, rates{rates}, STARTING_TYPE{starting_type} {
            for (int i=0; i<number_types; i++) {
                std::vector<int> cur_id_to_id_anc{};
                std::vector<short> cur_id_to_ancestor_type{};
                std::vector<int> cur_id_to_founder{};

                cur_id_to_id_anc.resize(MAX_TREE_SIZE_BY_TYPE[i]);
                cur_id_to_ancestor_type.resize(MAX_TREE_SIZE_BY_TYPE[i]);
                id_to_id_ancestry.push_back(cur_id_to_id_anc);
                id_to_ancestor_type.push_back(cur_id_to_ancestor_type);

                cur_id_to_founder.resize(MAX_TREE_SIZE_BY_TYPE[i]);
                id_to_founder.push_back(cur_id_to_founder);

                std::vector<double> cur_id_to_birth_time{};
                std::vector<double> cur_id_to_death_time{};
                cur_id_to_birth_time.resize(MAX_TREE_SIZE_BY_TYPE[i]);
                cur_id_to_death_time.resize(MAX_TREE_SIZE_BY_TYPE[i]);
                id_to_birth_time.push_back(cur_id_to_birth_time);
                id_to_death_time.push_back(cur_id_to_death_time);


                alive_pointer_left.push_back(0);
                alive_pointer_right.push_back(0);
            }
            id_to_founder[STARTING_TYPE][0]=0;
            add_to_alive(starting_type,starting_type,0);
        }

        Tree(const std::vector<int>& MAX_TREE_SIZE_BY_TYPE, const short STARTING_TYPE, std::vector<std::vector<int>>& id_to_id_ancestry, std::vector<std::vector<int>>& id_to_founder, std::vector<std::vector<short>>& id_to_ancestor_type, std::vector<std::vector<double>>& id_to_birth_time, std::vector<std::vector<double>>& id_to_death_time, const std::vector<std::vector<OffspringPossibility>>& offspring_distribution, std::vector<double>& rates, std::vector<int>& alive_pointer_left, std::vector<int>& alive_pointer_right): MAX_TREE_SIZE_BY_TYPE{MAX_TREE_SIZE_BY_TYPE}, STARTING_TYPE{STARTING_TYPE}, id_to_id_ancestry{id_to_id_ancestry}, id_to_founder{id_to_founder}, id_to_ancestor_type{id_to_ancestor_type}, id_to_birth_time{id_to_birth_time}, id_to_death_time{id_to_death_time}, offspring_distribution{offspring_distribution}, rates{rates}, alive_pointer_left{alive_pointer_left}, alive_pointer_right{alive_pointer_right} {
            number_types=MAX_TREE_SIZE_BY_TYPE.size();
        }

        void remove_from_alive(short type, int id) {
            int place_to_move_to{alive_pointer_left[type]};
            int place_to_move_to_value{id_to_id_ancestry[type][place_to_move_to]};
            int place_to_move_to_founder_value{id_to_founder[type][place_to_move_to]};

            id_to_id_ancestry[type][place_to_move_to]=id_to_id_ancestry[type][id];
            id_to_id_ancestry[type][id]=place_to_move_to_value;

            id_to_founder[type][place_to_move_to]=id_to_founder[type][id];
            id_to_founder[type][id]=place_to_move_to_founder_value;

            alive_pointer_left[type]++;

            int place_to_move_to_value_type{id_to_ancestor_type[type][place_to_move_to]};

            id_to_ancestor_type[type][place_to_move_to]=id_to_ancestor_type[type][id];
            id_to_ancestor_type[type][id]=place_to_move_to_value_type;
        }

        void add_to_alive(short cur_type, short ancestor_type, int ancestor_id) {
            int to_add_to{alive_pointer_right[cur_type]};
            id_to_id_ancestry[cur_type][to_add_to]=ancestor_id;
            id_to_ancestor_type[cur_type][to_add_to]=ancestor_type;

            if (id_to_ancestor_type[ancestor_type][ancestor_id]!=ancestor_type) id_to_founder[cur_type][to_add_to]=ancestor_id; //ancestor is a founder
            else id_to_founder[cur_type][to_add_to]=id_to_founder[ancestor_type][ancestor_id]; //ancestor is not a founder
            alive_pointer_right[cur_type]++;
            if (alive_pointer_right[cur_type]>=MAX_TREE_SIZE_BY_TYPE[cur_type]) std::cout << "ABOVE MAX SIZE: TYPE "+std::to_string(cur_type) << " | SIZE: " << alive_pointer_right[cur_type] << " | POP SIZE: " << alive_pointer_right[cur_type]-alive_pointer_left[cur_type];
        }

        template <typename T>
        void perform_save(Saver<T>& saver, short parent_type, double time, int event) {
            saver.add_save(parent_type,alive_pointer_left[parent_type]-1,time,offspring_distribution[parent_type][event].type_to_number,alive_pointer_left,alive_pointer_right);
        }

        void get_generations_at_saves(std::vector<SaveSnapshot> const& saves) {
            std::vector<std::vector<std::vector<std::vector<int>>>> output{};
            for (SaveSnapshot s: saves) {
                std::vector<std::vector<std::vector<int>>> current_gen_save{};
                for (int i=0;i<number_types;i++) {
                    std::vector<std::vector<int>> current_gen_save_current_type{};
                    current_gen_save_current_type.resize(s.alive_pointer_right[i]-s.alive_pointer_left[i]);
                    for (int j=0; j<s.alive_pointer_right[i]-s.alive_pointer_left[i];j++) {
                        std::vector<int> current_gen_save_current_type_current_cell{};
                        current_gen_save_current_type_current_cell.resize(number_types);
                        current_gen_save_current_type[j]=current_gen_save_current_type_current_cell;
                    }
                    current_gen_save.push_back(current_gen_save_current_type);

                }
                output.push_back(current_gen_save);
            }
        }



        
        void replace_cell_random_offspring(int type, int event, std::mt19937& rng) {
            std::uniform_int_distribution<> dist(alive_pointer_left[type],alive_pointer_right[type]-1);
            
            remove_from_alive(type,dist(rng));
            int i{0};
            for (int numb_off: offspring_distribution[type][event].type_to_number) {
                for (int j=0;j<numb_off;j++){
                    add_to_alive(i,type,alive_pointer_left[type]-1);
                }
                i++;
            }
        }

        template <typename T>
        std::unordered_map<std::string,IntDoubleVector> get_populations_at_last_save(Saver<T> const& saver, std::unordered_map<std::string,IntDoubleVector>& output) {

            auto last_save = saver.get_last_used_save();
            auto last_save_location = saver.get_last_used_save_location();

            if (!output.contains(saver.get_identifier_name())) output[saver.get_identifier_name()]=std::vector<T>{};
            for (short cur_type=0; cur_type<number_types; cur_type++) {
                if (!output.contains("population"+std::to_string(cur_type))) output["population"+std::to_string(cur_type)]=std::vector<int>{};
            }

            int cur_save{0};
            std::get<std::vector<T>>(output[saver.get_identifier_name()]).push_back(last_save_location);
            for (short cur_type=0; cur_type<number_types; cur_type++) {
                 std::get<std::vector<int>>(output["population"+std::to_string(cur_type)]).push_back(alive_pointer_right[cur_type]-alive_pointer_left[cur_type]);
            }
            return output;
        }

        template<typename T>
        std::unordered_map<std::string,IntDoubleVector> get_driver_SFS_at_last_save(Saver<T> const& saver, short type, std::unordered_map<std::string,IntDoubleVector>& output) {
            // Go from alive cells to founder
            //std::cout << "In get_driver_SFS_at_saves\n";
            if (!output.contains("type_"+std::to_string(type)+"_clone_size")) output["type_"+std::to_string(type)+"_clone_size"]=std::vector<int>();
            if (!output.contains("frequency")) output["frequency"]=std::vector<int>();
            if (!output.contains(saver.get_identifier_name())) output[saver.get_identifier_name()]=std::vector<T>();
            auto last_save = saver.get_last_used_save();
            auto last_save_location = saver.get_last_used_save_location();

            std::unordered_map<int,int> founder_of_type_to_number_descendents{};
            for (int cur_cell_id=alive_pointer_left[type]; cur_cell_id<alive_pointer_right[type]; cur_cell_id++) {
                //std::cout << (static_cast<double>(cur_cell_id-s.alive_pointer_left[type]))/(s.alive_pointer_right[type]-s.alive_pointer_left[type]) << '\n';
                /*int traverse_id=cur_cell_id;
                while (id_to_ancestor_type[type][traverse_id]==type && traverse_id!=0) {
                    traverse_id=id_to_id_ancestry[type][traverse_id];
                }
                if (!founder_of_type_to_number_descendents.contains(traverse_id)) founder_of_type_to_number_descendents[traverse_id]=0;
                founder_of_type_to_number_descendents[traverse_id]++;*/
                if (id_to_ancestor_type[type][cur_cell_id]!=type) {
                    founder_of_type_to_number_descendents[cur_cell_id]=1;
                } else {
                    int founder_id{id_to_founder[type][cur_cell_id]};
                    if (!founder_of_type_to_number_descendents.contains(id_to_founder[type][cur_cell_id])) founder_of_type_to_number_descendents[id_to_founder[type][cur_cell_id]]=0;
                    founder_of_type_to_number_descendents[id_to_founder[type][cur_cell_id]]++;
                }
            }

            std::unordered_map<int,int> size_of_clone_to_number{};
            for (auto const& [_,value] : founder_of_type_to_number_descendents) {
                if (!size_of_clone_to_number.contains(value)) size_of_clone_to_number[value]=0;
                size_of_clone_to_number[value]++;
            }

            for (auto const & [k,v] : size_of_clone_to_number) {
                std::get<std::vector<int>>(output["type_"+std::to_string(type)+"_clone_size"]).push_back(k);
                std::get<std::vector<int>>(output["frequency"]).push_back(v);
                std::get<std::vector<T>>(output[saver.get_identifier_name()]).push_back(last_save_location);
            }
            return output;
        }



        /*void get_clonal_ancestry_at_saves(std::vector<std::unordered_map<int,int>>& id_to_id_founders,std::vector<std::unordered_map<int,short>>& id_to_prev_type_founders, std::vector<std::unordered_map<int,double>>& id_to_initiation_time) {
            id_to_id_founders[STARTING_TYPE][0]=0;
            id_to_prev_type_founders[STARTING_TYPE][0]=STARTING_TYPE;
            id_to_initiation_time[STARTING_TYPE][0]=0;
            for (SaveSnapshot s: saves) {
                id_to_initiation_time[s.type_of_cell_of_interest][s.id_of_cell_of_interest]=s.time;
                for (int type=0; type<number_types; type++) {
                    for (int cur_cell=s.alive_pointer_left[type]; cur_cell<s.alive_pointer_right[type]; cur_cell++) {
                        int cell_traversal{cur_cell};
                        int founder_cell{-1};
                        int founder_type_next{type};
                        while (!id_to_id_founders[founder_type_next].contains(cell_traversal)) {
                            if (id_to_ancestor_type[founder_type_next][cell_traversal]!=founder_type_next) {
                                if (founder_cell!=-1) {
                                    id_to_id_founders[founder_type_next][founder_cell]=cell_traversal;
                                    id_to_prev_type_founders[founder_type_next][founder_cell]=id_to_ancestor_type[founder_type_next][cell_traversal];
                                }
                                founder_cell=cell_traversal;
                                cell_traversal=id_to_id_ancestry[founder_type_next][cell_traversal];
                                founder_type_next=id_to_ancestor_type[founder_type_next][cell_traversal];
                            }
                            else {
                                cell_traversal=id_to_id_ancestry[founder_type_next][cell_traversal];
                            }
                            
                        }
                        if (founder_cell!=-1) {
                            id_to_id_founders[founder_type_next][founder_cell]=cell_traversal;
                            id_to_prev_type_founders[founder_type_next][founder_cell]=id_to_ancestor_type[founder_type_next][cell_traversal];
                        }

                    }
                }
            }
        }*/

        int get_number_types() {
            return number_types;
        }
        int get_population(short type) {
            return alive_pointer_right[type] - alive_pointer_left[type];
        }
        int get_tree_size(short type) {
            return alive_pointer_right[type];
        }
        int get_max_tree(short type) {
            return MAX_TREE_SIZE_BY_TYPE[type];
        }
        double get_rate(short type) {
            return rates[type];
        }
        const std::vector<OffspringPossibility>& get_offspring_distribution(short type) {
            return offspring_distribution[type];
        }
        const std::vector<std::vector<OffspringPossibility>>& get_offspring_distribution() {
            return offspring_distribution;
        }

        Tree obtain_tree_from_saved_old_to_new_offset(const std::vector<std::map<int,int>>& old_id_to_new_id_offset, const std::vector<int>& offsets) {

            std::vector<int> new_MAX_TREE_SIZE_BY_TYPE{};
            std::vector<std::vector<int>> new_id_to_id_ancestry{};
            std::vector<std::vector<int>> new_id_to_id_founder{};
            std::vector<int> new_alive_pointer_left{};
            std::vector<int> new_alive_pointer_right{};
            std::vector<std::vector<short>> new_id_to_ancestor_type{};
            std::vector<std::vector<double>> new_id_to_birth_time{};
            std::vector<std::vector<double>> new_id_to_death_time{};

            for (int i=0; i<number_types; i++) {
                new_MAX_TREE_SIZE_BY_TYPE.push_back(old_id_to_new_id_offset[i].size());
                std::vector<int> cur_new_id_to_id_ancestry{};
                std::vector<int> cur_new_id_to_id_founder{};
                std::vector<short> cur_new_id_to_ancestor_type{};
                std::vector<double> cur_new_id_to_birth_time{};
                std::vector<double> cur_new_id_to_death_time{};

                cur_new_id_to_id_ancestry.resize(new_MAX_TREE_SIZE_BY_TYPE[i]);
                cur_new_id_to_id_founder.resize(new_MAX_TREE_SIZE_BY_TYPE[i]);
                cur_new_id_to_ancestor_type.resize(new_MAX_TREE_SIZE_BY_TYPE[i]);
                cur_new_id_to_birth_time.resize(new_MAX_TREE_SIZE_BY_TYPE[i]);
                cur_new_id_to_death_time.resize(new_MAX_TREE_SIZE_BY_TYPE[i]);

                new_id_to_id_ancestry.push_back(cur_new_id_to_id_ancestry);
                new_id_to_id_founder.push_back(cur_new_id_to_id_founder);
                new_id_to_ancestor_type.push_back(cur_new_id_to_ancestor_type);
                new_id_to_birth_time.push_back(cur_new_id_to_birth_time);
                new_id_to_death_time.push_back(cur_new_id_to_death_time);

                new_alive_pointer_right.push_back(new_MAX_TREE_SIZE_BY_TYPE[i]);
                new_alive_pointer_left.push_back(offsets[i]);
            }

            for (int i=0; i<number_types; i++) {
                int offset{offsets[i]};
                for (auto const& [k,v]: old_id_to_new_id_offset[i]) {
                    int anc_type{id_to_ancestor_type[i][k]};
                    int anc_id{id_to_id_ancestry[i][k]};

                    new_id_to_id_ancestry[i][v+offset]=old_id_to_new_id_offset[anc_type].at(anc_id)+offset;
                    new_id_to_id_founder[i][v+offset]=id_to_founder[i][k];
                    new_id_to_ancestor_type[i][v+offset]=anc_type;
                    new_id_to_birth_time[i][v+offset]=id_to_birth_time[i][k];
                    new_id_to_death_time[i][v+offset]=id_to_death_time[i][k];
                }
            }

            Tree out{new_MAX_TREE_SIZE_BY_TYPE,STARTING_TYPE,new_id_to_id_ancestry,new_id_to_id_founder,new_id_to_ancestor_type,new_id_to_birth_time,new_id_to_death_time,offspring_distribution,rates,new_alive_pointer_left,new_alive_pointer_right};
            return out;

        }

        Tree sample_subtree_binomial(const std::vector<double>& sampling_probabilities_by_type, std::mt19937& rng) {
            std::vector<int> reserved_size_alive{};
            std::vector<std::vector<int>> id_sample_alive{};
            std::vector<std::bernoulli_distribution> bernoullis{};
            std::vector<std::map<int,int>> old_id_to_new_id_offset{};
            //Initialization
            for (int i=0; i<sampling_probabilities_by_type.size(); i++) {
                reserved_size_alive.push_back(static_cast<int>((alive_pointer_right[i]-alive_pointer_left[i]) * sampling_probabilities_by_type[i]+2.5*sqrt(alive_pointer_right[i]-alive_pointer_left[i])));
                std::vector<int> cur_type_sample_id_alive{};
                cur_type_sample_id_alive.reserve(reserved_size_alive[i]);
                id_sample_alive.push_back(cur_type_sample_id_alive);
                bernoullis.push_back(std::bernoulli_distribution{sampling_probabilities_by_type[i]});
                old_id_to_new_id_offset.push_back(std::map<int,int>{});
            }

            for (int i=0; i<number_types; i++) {
                int new_ind_off{0};
                for (int cur_ind=alive_pointer_left[i]; cur_ind < alive_pointer_right[i]; cur_ind++) {
                    if (bernoullis[i](rng)) {
                        id_sample_alive[i].push_back(cur_ind);
                        old_id_to_new_id_offset[i][cur_ind]=new_ind_off;
                        new_ind_off++;
                    }
                }
            }

            
            std::vector<int> offsets{};
            int last_old_id_of_starting_type_seen{};
            for (int i=0; i<number_types; i++) {
                int cur_new_id{-1};
                int cur_number_cells_in_anc{0};
                for (int j=0; j<id_sample_alive[i].size(); j++) {
                    int cur_id{id_sample_alive[i][j]};
                    int cur_type{i};
                    old_id_to_new_id_offset[i][cur_id]=j;
                    while (id_to_id_ancestry[cur_type][cur_id]!=cur_id || id_to_ancestor_type[cur_type][cur_id]!=cur_type) {
                        cur_id=id_to_id_ancestry[cur_type][cur_id];
                        cur_type=id_to_ancestor_type[cur_type][cur_id];
                        if (old_id_to_new_id_offset[cur_type].contains(cur_id)) break;
                        old_id_to_new_id_offset[cur_type][cur_id]=cur_new_id;
                        if (i==STARTING_TYPE) last_old_id_of_starting_type_seen=cur_id;
                        cur_new_id--;
                        cur_number_cells_in_anc++;
                    }
                }
                offsets.push_back(-cur_new_id-1);
            }

            //Make the root the earliest type.
            old_id_to_new_id_offset[STARTING_TYPE][last_old_id_of_starting_type_seen]=old_id_to_new_id_offset[STARTING_TYPE][0];
            old_id_to_new_id_offset[STARTING_TYPE][0]=-offsets[STARTING_TYPE];

            return obtain_tree_from_saved_old_to_new_offset(old_id_to_new_id_offset,offsets);
        }

        Tree binomial_rarify(std::vector<std::vector<double>> rarify_edge_probabilities_by_type, std::mt19937& rng) {
            std::vector<std::map<int,bool>> keep_founders_by_type{};

            std::vector<std::vector<std::bernoulli_distribution>> bernoullis{};
            std::vector<std::map<int,int>> old_id_to_new_id_offset{};

            //Initialize Bernoullis
            for (int i=0; i<number_types; i++) {
                bernoullis.push_back(std::vector<std::bernoulli_distribution>{});
                keep_founders_by_type.push_back(std::map<int,bool>{});
                for (int j=0; j<number_types; j++) {
                    bernoullis[i].push_back(std::bernoulli_distribution(rarify_edge_probabilities_by_type[i][j]));
                }
                old_id_to_new_id_offset.push_back(std::map<int,int>{});
            }
            keep_founders_by_type[STARTING_TYPE][0]=true;

            for (int i=0; i<number_types; i++) {
                for (int j=alive_pointer_right[i]-1; j>=0; j--) {
                    int cur_founder_id{j};
                    int cur_founder_type{i};
                    int first_founder_id_false_in_line{-1};
                    int first_founder_type_false_in_line{-1};

                    if (id_to_ancestor_type[i][j]==i) {
                        cur_founder_id=id_to_founder[i][j];
                    }

                    while (!keep_founders_by_type[cur_founder_type].contains(cur_founder_id)) {
                        keep_founders_by_type[cur_founder_type][cur_founder_id]=bernoullis[id_to_ancestor_type[cur_founder_type][cur_founder_id]][cur_founder_type](rng);
                        if (!keep_founders_by_type[cur_founder_type][cur_founder_id]) {
                            first_founder_id_false_in_line=cur_founder_id;
                            first_founder_type_false_in_line=cur_founder_type;
                        }

                        int temp_founder_id{cur_founder_id};
                        cur_founder_id=id_to_founder[cur_founder_type][temp_founder_id];
                        cur_founder_type=id_to_ancestor_type[cur_founder_type][temp_founder_id];
                    }
                    if (!keep_founders_by_type[cur_founder_type][cur_founder_id]) {
                        first_founder_id_false_in_line=cur_founder_id;
                        first_founder_type_false_in_line=cur_founder_type;
                    }

                    cur_founder_id=j;
                    cur_founder_type=i;
                    if (id_to_ancestor_type[i][j]==i) {
                        cur_founder_id=id_to_founder[i][j];
                    }
                    while (cur_founder_id!=first_founder_id_false_in_line || cur_founder_type!=first_founder_type_false_in_line) {
                        keep_founders_by_type[cur_founder_type][cur_founder_id]=false;
                        int temp_founder_id{cur_founder_id};
                        cur_founder_id=id_to_founder[cur_founder_type][temp_founder_id];
                        cur_founder_type=id_to_ancestor_type[cur_founder_type][temp_founder_id];
                    }
                }
            }
            //With kept founders, iterate and keep necessary cells.
            std::vector<int> offsets{};
            for (int i=0; i<number_types; i++) {
                int cur_alive_cell_ind_offset{0};
                int cur_anc_cell_ind_offset{-1};
                for (int cur_cell_id=alive_pointer_right[i]-1; cur_cell_id>=0; cur_cell_id--) {
                    int cur_founder{cur_cell_id};
                    if (id_to_ancestor_type[i][cur_cell_id]==i) cur_founder=id_to_founder[i][cur_cell_id];

                    if (!keep_founders_by_type[i][cur_founder]) continue;

                    if (cur_cell_id>=alive_pointer_left[i]) {
                        old_id_to_new_id_offset[i][cur_cell_id]=cur_alive_cell_ind_offset;
                        cur_alive_cell_ind_offset++;
                    } else {
                        old_id_to_new_id_offset[i][cur_cell_id]=cur_anc_cell_ind_offset;
                        cur_anc_cell_ind_offset--;
                    }
                }
                offsets.push_back(-cur_anc_cell_ind_offset-1);
            }
            return obtain_tree_from_saved_old_to_new_offset(old_id_to_new_id_offset,offsets);


        }
};

template <typename T>
class Saver {
    protected:
        const std::vector<T> SAVE_LOCATIONS{};
        std::vector<SaveSnapshot> saves;
        int current_save_pointer{0};
        int last_used_save_pointer{0};
    public:
        Saver(std::vector<T> save_locations): SAVE_LOCATIONS{save_locations}{}

        void move_to_next_save() {
            current_save_pointer++;
        }

        bool check_last_used_save() {
            if (last_used_save_pointer>=current_save_pointer) return false;
            return true;
        }

        void increment_last_used_save() {
            last_used_save_pointer++;
        }

        virtual void add_save(short parent_type, int parent_id, double time, std::vector<int> const& number_of_offspring_of_type, std::vector<int> const& alive_pointer_left, std::vector<int> const& alive_pointer_right)=0;
        std::vector<T> get_save_locations(){return SAVE_LOCATIONS;}
        virtual std::string get_identifier_name()const =0;

        std::vector<SaveSnapshot> const& get_saves() const {
            return saves;
        }

        SaveSnapshot const& get_last_used_save() const {
            return saves[last_used_save_pointer];
        }

        int get_last_used_save_location_index() const {
            return last_used_save_pointer;
        }

        T const& get_last_used_save_location() const {
            return SAVE_LOCATIONS[last_used_save_pointer];
        }

        std::vector<T> const& get_save_locations() const {
            return SAVE_LOCATIONS;
        }
};

class SizeSaver: public Saver<int> {
    public:
        inline static const std::string IDENTIFIER_NAME="population";
        SizeSaver(std::vector<int> save_locations): Saver(save_locations){}

        void add_save(short parent_type,int parent_id, double time, std::vector<int> const& number_of_offspring_of_type, std::vector<int> const& alive_pointer_left, std::vector<int> const& alive_pointer_right) override {
            if (current_save_pointer<0) {
                return;
            }

            int sum_total_sizes{0};
            for (int i=0; i<alive_pointer_left.size(); i++) {
                sum_total_sizes+=alive_pointer_right[i]-alive_pointer_left[i];
            }
            if (current_save_pointer>=SAVE_LOCATIONS.size() || sum_total_sizes<SAVE_LOCATIONS[current_save_pointer]) return;
            std::vector<int> total_sizes{alive_pointer_right};
            for (int i=0; i<alive_pointer_left.size(); i++){
                total_sizes[i]-=alive_pointer_left[i];
            }

            while (current_save_pointer<SAVE_LOCATIONS.size() && sum_total_sizes>=SAVE_LOCATIONS[current_save_pointer]) {
                SaveSnapshot cur_save{alive_pointer_left,alive_pointer_right,time,parent_type,parent_id};
                saves.push_back(cur_save);
                current_save_pointer++;
            }
        }

        std::string get_identifier_name() const override {
            return SizeSaver::IDENTIFIER_NAME;
        }


};

class TreeSizeSaver: public Saver<int> {
    public:
        inline static const std::string IDENTIFIER_NAME="tree size";
        TreeSizeSaver(std::vector<int> save_locations): Saver(save_locations){}

        void add_save(short parent_type,int parent_id, double time, std::vector<int> const& number_of_offspring_of_type, std::vector<int> const& alive_pointer_left, std::vector<int> const& alive_pointer_right) override {
            if (current_save_pointer<0) {
                return;
            }

            int sum_total_sizes{0};
            for (int i=0; i<alive_pointer_left.size(); i++) {
                sum_total_sizes+=alive_pointer_right[i];
            }
            if (current_save_pointer>=SAVE_LOCATIONS.size() || sum_total_sizes<SAVE_LOCATIONS[current_save_pointer]) return;
            std::vector<int> total_sizes{alive_pointer_right};

            while (current_save_pointer<SAVE_LOCATIONS.size() && sum_total_sizes>=SAVE_LOCATIONS[current_save_pointer]) {
                SaveSnapshot cur_save{alive_pointer_left,alive_pointer_right,time,parent_type,parent_id};
                saves.push_back(cur_save);
                current_save_pointer++;
            }
        }

        std::string get_identifier_name() const override {
            return TreeSizeSaver::IDENTIFIER_NAME;
        }


};

template<int TYPE>
class TypeSizeSaver: public Saver<int> {
    public:
        inline static const std::string IDENTIFIER_NAME="population-"+std::to_string(TYPE);
        TypeSizeSaver(std::vector<int> save_locations): Saver(save_locations){}

        void add_save(short parent_type,int parent_id, double time, std::vector<int> const& number_of_offspring_of_type, std::vector<int> const& alive_pointer_left, std::vector<int> const& alive_pointer_right) override {
            if (current_save_pointer<0) {
                return;
            }

            int sum_total_sizes{alive_pointer_right[TYPE]-alive_pointer_left[TYPE]};
            if (current_save_pointer>=SAVE_LOCATIONS.size() || sum_total_sizes<SAVE_LOCATIONS[current_save_pointer]) return;
            std::vector<int> total_sizes{alive_pointer_right};

            while (current_save_pointer<SAVE_LOCATIONS.size() && sum_total_sizes>=SAVE_LOCATIONS[current_save_pointer]) {
                SaveSnapshot cur_save{alive_pointer_left,alive_pointer_right,time,parent_type,parent_id};
                saves.push_back(cur_save);
                current_save_pointer++;
            }
        }

        std::string get_identifier_name() const override {
            return TypeSizeSaver::IDENTIFIER_NAME;
        }
};

class TimeSaver: public Saver<double> {
    public:
        inline static const std::string IDENTIFIER_NAME="time";
        TimeSaver(std::vector<double> save_locations): Saver(save_locations){}

        void add_save(short parent_type,int parent_id, double time, std::vector<int> const& number_of_offspring_of_type, std::vector<int> const& alive_pointer_left, std::vector<int> const& alive_pointer_right) override {
            if (current_save_pointer<0) {
                return;
            }

            while (current_save_pointer<SAVE_LOCATIONS.size() && time>=SAVE_LOCATIONS[current_save_pointer]) {
                std::vector<int> prev_alive_pointer_right{alive_pointer_right};
                std::vector<int> prev_alive_pointer_left{alive_pointer_left};
                for (short cur_type=0; cur_type<alive_pointer_right.size(); cur_type++){
                    prev_alive_pointer_right[cur_type]-=number_of_offspring_of_type[cur_type];
                    if (cur_type==parent_type) prev_alive_pointer_left[cur_type]-=1;
                }
                //std::cout<<prev_alive_pointer_right[0]-alive_pointer_left[0]<<'\n';
                SaveSnapshot cur_save{prev_alive_pointer_left,prev_alive_pointer_right,SAVE_LOCATIONS[current_save_pointer],parent_type, parent_id};
                saves.push_back(cur_save);
                current_save_pointer++;

            }
        }

        std::string get_identifier_name() const override {
            return TimeSaver::IDENTIFIER_NAME;
        }
};

class MutantOffspringSaver: public Saver<std::monostate> {
    public:
        MutantOffspringSaver(): Saver({}){}

        void add_save(short parent_type,int parent_id, double time, std::vector<int> const& number_of_offspring_of_type, std::vector<int> const& alive_pointer_left, std::vector<int> const& alive_pointer_right) override {
            for (short cur_type=0; cur_type<number_of_offspring_of_type.size(); cur_type++) {
                if (parent_type==cur_type) continue;
                for (int j=0; j<number_of_offspring_of_type[cur_type]; j++) {
                    int cur_id=alive_pointer_right[cur_type]-j-1;
                    SaveSnapshot cur_save{alive_pointer_left,alive_pointer_right,time,cur_type,cur_id};
                    saves.push_back(cur_save);
                }
            }
        }

        std::string get_identifier_name() const override {
            return "";
        }
};



#endif