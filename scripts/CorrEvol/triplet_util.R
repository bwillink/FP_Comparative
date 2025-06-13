# probabilities of next event type
make_P = function(Q) {
    P = Q
    for (i in 1:nrow(Q)) {
        P[i,] = Q[i,] / -Q[i,i]
        P[i,i] = 0
    }
    return(P)
}

# connectivity
make_Q = function(params, graph) {
    
    ns = 4
    nt = length(graph)
    
   # raw biome shift rates
    br = matrix(0, nrow=ns, ncol=ns)
    br[1,2] = params$rate_gain_A_when_B0
    br[1,3] = params$rate_gain_B_when_A0
    br[2,4] = params$rate_gain_B_when_A1
    br[3,4] = params$rate_gain_A_when_B1
    br[2,1] = params$rate_loss_A_when_B0
    br[3,1] = params$rate_loss_B_when_A0
    br[4,2] = params$rate_loss_B_when_A1
    br[4,3] = params$rate_loss_A_when_B1
    
    Q = matrix(NA, nrow=ns, ncol=ns)
    #rownames(Q) <- colnames(Q) <- abnames
    
    for (bi in 1:ns) {
      for (bj in 1:ns) {
        biome_diff = F
        
        if (bi == bj) {
          # do nothing
        } else {
          biome_diff = T
        }
        
        g_biome = graph
        
        Q[bi, bj] = 0
        if (biome_diff) {
          # biome shift
          Q[bi, bj] = br[bi, bj]
        }
      }
    }
    
    for (i in 1:ns) {
        Q[i,i] = -sum(Q[i,])
    }
    
    return(Q)
}

make_triplets = function(x, params, graph, subroot=F) {
    
    # get root node and state
    idx_root = max(x$node_index)
    x_root = x[ x$parent_index==idx_root, ]$start_state[1]
    
    # get root state Q matrix
    Q_root = make_Q( params, graph )
    
    # turn into probability matrix for next event being type ij
    P_root = make_P( Q_root )
    
    # get root state probabilities
    rf_lbl = paste("rf.",1:ncol(Q_root),".",sep="")
    rf_prob = params[,c(rf_lbl)]
    
    # compute the probability of the previous state being X_subroot=i given you end in state X_root=j
    
    # Define:
    #   p1 : P( X_subroot=i | X_root=j, Q(m_root) ) -- want, but don't have
    #   p2 : P( X_root=j | X_subroot=i, Q(m_root) ) -- prob event of j | i
    #   p3 : pi( X_subroot=i | Q(m_root) ) -- stationary prob i
    #   p4 : pi( X_root=j | Q(m_root) -- stationary prob j
    #   p1 = p2 * p3 / p4 -- Bayes
    
    # fill in our probabilities for starting in i given we end in j
    rev_prob = rep(0, 4)
    for (i in 1:length(rf_prob)) {
        rev_prob[i] = P_root[i, x_root] * rf_prob[i] / rf_prob[x_root]
    }
    rev_prob = unlist(rev_prob)
    rev_prob = rev_prob / sum(rev_prob)
    
    # sample subroot state
    x_subroot = -1
    if (subroot) {
        x_subroot = sample(1:4, size=1, prob=rev_prob, replace=TRUE)
    }
    state_root = c(idx_root, -1, x_subroot, x_root)
    count_root = 0
    
    # start recursion
    events = make_triplet_recursion(x, idx_root, state_root, count_root)
    
    # return count mtx
    return(events)
}

make_triplet_df = function(d,params,graph,subroot=F) {
    # loop over iterations
    iterations = sort(unique(d$iteration))
    n_it = length(iterations)
    
    evt = NULL
    df_list = list()
    df_list_idx = 1
    for ( i in 1:n_it ) {
        # get iteration
        it = iterations[i]
        # get data for iteration
        dtmp = d[d$iteration==it, ]
        # get parameters for iteration
        partmp = params[ params$Iteration==it, ]
        # make event triplets for iteration
        evttmp = data.frame(it, make_triplets( dtmp, partmp, graph, subroot ))
        # add event triplets to data frame
        #evt = rbind(evt, evttmp)
        df_list[[df_list_idx]] = evttmp
        df_list_idx = df_list_idx + 1
    }
    evt = rbindlist(df_list)
    
    # return data frame
    df = data.frame(evt); colnames(df) = c("iteration","node_index","S1","S2","S3","age"); rownames(df)=NULL
    return(df)
}

make_dat_triplet_raw = function(fn, n_states, f_burn=0.0, thinby=1) { 
    
    # files
    #out_fp = paste(fp, "output_trim/", sep="")
    #fn = paste(out_fp, base_fn, ".history.tsv", sep="")
    
    # read in data file
    stoch = read.csv(fn, sep="\t", stringsAsFactors=F)
    stoch = stoch[ stoch$transition_type != "cladogenetic", ]
    stoch$transition_time[ stoch$transition_type=="no_change" ] = stoch$branch_start_time[ stoch$transition_type=="no_change" ]
    
    # filter MCMC samples
    iterations = unique(stoch$iteration)
    n_burn = max(1, f_burn*length(iterations))
    iterations = iterations[n_burn:length(iterations)]
    iterations = iterations[ seq(1, length(iterations), by=thinby)  ]
        
    # get branch indexing
    branches = 1:max(unique(stoch$parent_index), na.rm=T)
    
    # STAGE 1, gather data
    # loop over iterations
    df_list = list()
    df_list_idx = 1
   
    for (i in 1:length(iterations)) {
        
        # get biome and biogeography stochastic mappings per iteration
        it = iterations[i]
        cat("Stage 1, processing iteration ",it," / ", max(iterations), "\n", sep="")
        sample = stoch[ stoch$iteration==it, ]
        
        # loop over branches
        br_list = list()
        br_list_idx = 1
        
        dtmp = data.frame(stringsAsFactors = F)
        for (j in 1:length(branches)) {
            
            # get biome and biogeography stochastic mappings per branch
            nd_idx = branches[j]
            branch = sample[ sample$node_index==nd_idx, ]
            
            # interleave biome and biogeography stochastic mappings
            branch_bg_biome = branch
            branch_bg_biome$start_state = branch_bg_biome$start_state + 1
            branch_bg_biome$end_state = branch_bg_biome$end_state + 1
            
            dtmp = as.data.frame(branch_bg_biome, stringsAsFactors=F)
            br_list[[br_list_idx]] = dtmp
            br_list_idx = br_list_idx + 1
            
        }
        df_br = rbindlist(br_list)
        df_list[[df_list_idx]] = df_br
        df_list_idx = df_list_idx + 1
        
    }
    stoch_bg_biome = rbindlist(df_list)
    
    return(stoch_bg_biome)
}

make_triplet_recursion = function(x, idx, anc_state, count) {
  events = matrix(NA, nrow=0, ncol=5)
  
  # get node
  x_nd = x[ x$node_index==idx, ]
  
  # collect events if something happened on branch
  nr = nrow(x_nd)
  
  if (nr > 0 && x_nd$transition_type[1] != "no_change" && !is.na(x_nd$parent_index[1]) ) {
    
    # loop over events
    for (i in 1:nr) {
      
      # get prev and curr state
      state_0 = anc_state[3]
      state_1 = x_nd$start_state[i]
      state_2 = x_nd$end_state[i]
      age = x_nd$transition_time[i]
      
      # update triplet vector
      anc_state = c(idx, state_0, state_1, state_2, age)
      
      # append to event dataframe
      events = rbind(events, anc_state)
    }
  }
  
  # recurse to children
  ch1_events = NULL
  ch1_idx = x_nd$child1_index[1]
  #print(ch1_idx)
  if (!is.na(ch1_idx)) {
    ch1_events = make_triplet_recursion(x, ch1_idx, anc_state, count)
  }
  ch2_events = NULL
  ch2_idx = x_nd$child2_index[1]
  #print(ch2_idx)
  if (!is.na(ch2_idx)) {
    ch2_events = make_triplet_recursion(x, ch2_idx, anc_state, count)
  }
  
  # collect events
  events = rbind(events, ch1_events, ch2_events)
  
  return(events)
}