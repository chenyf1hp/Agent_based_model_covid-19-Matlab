    for sym_count = 1:length(sym_index)
        infectious_flag = 3;
        origin_dis = people(sym_index(sym_count),10);
        sym_contact_n = N_contact(infectious_count);
        source_uid = people(sym_index(sym_count),1);
        dist_u = rand(sym_contact_n,1);
        [contact_index, contact_uid, infected_index, infected_uid]  = get_contact_uid(people,cum_p_contact,dis_index_pool,dist_u,contact_n,origin_dis,infectious_flag,r,beta);
        sym_source = [infected_uid,repmat(source_uid,length(infected_uid),1)];
        transmission_source = [transmission_source;sym_source];
        if isempty(transmission_chain)
            transmission_chain(infectious_count).source_uid = source_uid;
            transmission_chain(infectious_count).contact_uid = contact_uid;
            transmission_chain(infectious_count).infected_uid = infected_uid;
        else
            exist_source_uid = cat(1,transmission_chain.source_uid);
            if find(exist_source_uid==source_uid)
                source_index = find(exist_source_uid==source_uid);
                chain_contact_uid = [transmission_chain(source_index).contact_uid;contact_uid];
                chain_contact_uid = unique(chain_contact_uid);
                transmission_chain(source_index).contact_uid = chain_contact_uid;
                chain_infected_uid = [transmission_chain(source_index).infected_uid;infected_uid];
                chain_infected_uid = unique(chain_infected_uid);
                transmission_chain(source_index).infected_uid = chain_infected_uid;
            else
                transmission_chain(infectious_count).source_uid = source_uid;
                transmission_chain(infectious_count).contact_uid = contact_uid;
                transmission_chain(infectious_count).infected_uid = infected_uid;                
            end
        end