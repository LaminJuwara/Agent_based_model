#######################################################################################
# Generating uniformly distributed simulated data within
# 10% of default parameter
# 
#######################################################################################

pa = 1; pad = 0.50;          # default parameter values
p1 = 0.5; p3= 0.5 ;       

# Used the default parameter in miao's paper for exponential times
mu_p0 = 40;  mu_d0 = 60; mu_p = 8; mu_d = 60; mu_a = 16;

P = c(pa, pad, p1, p3, mu_a, mu_p0, mu_d0, mu_p, mu_d ) # required parameters
P_rand <- runif(length(P), min = 0.9*P , max = 1.1*P)
t_0 = 0 ; 
N_init = 50000;  # unactivated cells

# creating containers for generations with activation a, death d, and proliferation p.
gen0d = 0  ; gen0p=0  ;   gen1d = 0  ; gen1p = 0  ; gen2d = 0 ; gen2p=0 ;
gen3d = 0  ; gen3p=0  ;   gen4d = 0  ; gen4p = 0  ; gen5d = 0 ; gen5p=0 ; 
gen6d = 0  ; gen6p=0  ;   gen7d = 0  ; gen7p = 0  ; gen8d = 0 ; gen8p=0 ;
gen9d = 0  ; gen9p=0  ;   gen10d = 0 ; gen10p = 0 ; gen11d = 0; gen11p=0;
gen12d = 0 ; gen12p=0 ;   gen13d = 0 ; gen13p = 0 ; gen14d = 0; gen14p=0;
gen15d = 0 ; gen15p=0 ;   gen16d = 0 ; gen16p = 0 ; gen17d = 0; gen17p=0;
gen18d = 0 ; gen18p=0 ;   gen19d = 0 ; gen19p = 0 ; gen20d = 0; gen20p=0;


# Initializing 


# construct a  vector r_0 of uniform random length N_int 
ra <- runif(N_init , min = 0, max = 1)  
g_aindex = which(ra <= P[1])       # Proportion of cells programmed for activation N_0*Pa
g_0a = ra[g_aindex]                # cells programmed to be activated are stored in the vector

gen0a <- length(g_0a)  # store number of cells programmed to be activated
N_0r <- N_init - gen0a # number of cells programmed to remain resting-- no use for us

# determine the activated cells programmed for proliferation
rp = runif(gen0a , min = 0, max = 1)
g_pindex <- which(rp < (1-P_rand[2]))
g_0p = rp[g_pindex]     
# determine cells programmed for dying
g_dindex <- which(rp >= (1-P_rand[2]))
g_0d <- rp[g_dindex]


gen0p <- length(g_0p)  # store number of cells programmed to proliferate .i.e. N_init*Pa*(1-Pad)
gen0d <- length(g_0d)  # store number of cells programmed to be die

# only one event takes place at a time.
# To determine the first event that will take place
#############################################################################################
# determine rates for cells that will proliferate with mean mu_p0
lip0 = rexp(gen0p , rate = 1/P_rand[6])

# determine rates for cells that will die with mean mu_d0
lid0 = rexp(gen0d , rate = 1/P_rand[7])

###########################################################################################
tau_0 = which.min(c(lip0, lid0))           # determine the least time tau_0 index
delta_t0 = min(c(lip0, lid0))              # least time

# need to determine the first event that occured (proliferation or death)
if (tau_0 <= length(lip0) ) {   # if the least time is from cells programmed to divide
  gen0p = gen0p - 1
  # 2 two daughter cells are produced 
  U = runif(1, min=0, max=1)
  if (U <= P_rand[4]^2) { # both daughter cell are programmed to proliferate
    gen1p = gen1p + 2
  } 
  else if(U >= (1-(1-P_rand[4])^2)) { # both daughter cells are  programmed to die
    gen1d = gen1d + 2 
  }
  else{ # one is programmed to proliferate and one is programmed to die with prob 2P[death]P[div]
    gen1p = gen1p + 1
    gen1d = gen1d + 1
  }
} else {  # least time delta_to is from cells programmed to die in the zeroth generation
    gen0d = gen0d -1   # reduce the number of cells programmed to die
}
t_0 = t_0 + delta_t0       # update the time

#############################################################################################
# Now we have atleast two cell in generation one programmed to proliferate or die
#############################################################################################

# determine where the next event will occur i.e proliferation or death in all generations

Tt = t_0                    # renaming time steps from the initial event
Tmax = 120                  # set the max time
t.pts <- seq(0,120,12)      # set the time points in [0,120]
q <- 2
tpts <-t.pts[q]
lresult <- c()              # list for results and Matrix for storing cell counts
Ncell <- matrix( rep(0, (length(t_pts)*(ngen+1))), nrow = length(t_pts), ncol = ngen+1)



i = 2  # index for later results

while (Tt < Tmax) {

  # determine rates for cells that will proliferate with mean mu_p0
  li.p0 = rexp(gen0p, rate = 1/P_rand[6])
  # determine rates for cells that will die with mean mu_d0
  li.d0 = rexp(gen0d , rate = 1/P_rand[7])

  # determine rates for cells that will proliferate with mean mu_p
  # determine rates for cells that will die with mean mu_d
  for(j in 1:20){
    assign(paste('li.p',j,sep=''), rexp(get(paste('gen',j,'p',sep='')), rate = 1/P_rand[8]))
    assign(paste('li.d',j,sep=''), rexp(get(paste('gen',j,'d',sep='')), rate = 1/P_rand[9]))
  }
  
  # create a list for the minimum time in each division class
  dlist <- list(min(li.p0), min(li.d0), min(li.p1), min(li.d1), min(li.p2), min(li.d2), min(li.p3), min(li.d3),
                     min(li.p4), min(li.d4), min(li.p5), min(li.d5), min(li.p6), min(li.d6), min(li.p7),
                     min(li.d7), min(li.p8), min(li.d8), min(li.p9), min(li.d9),min(li.p10), min(li.d10), 
                     min(li.p11), min(li.d11), min(li.p12), min(li.d12), min(li.p13),
                     min(li.d13), min(li.p14), min(li.d14),min(li.p15), min(li.d15), min(li.p16), min(li.d16),
                     min(li.p17), min(li.d17), min(li.p18), min(li.d18),
                     min(li.p19), min(li.d19), min(li.p20), min(li.d20))
  clist <- c()
  for (i in 1:42) {  # create a container for the minimum times in each division class
    clist <- c(clist, dlist[[i]])
  }
  delta_t <- min(clist)
  x = which.min(clist) 
  
  
  
  if (x<=40) {
    if ( (x%% 2) == 0 ) { # check if the minimum index is even or odd to determine position
      k = (0.5*x) - 1
    } else{               # do this for odd position 
      k = (0.5*x) - 0.5   
    }

    # we determine the generation and specific death/proliferating class 
    # we update number of cells accordingly

    if (delta_t %in% get(paste('li.p',k,sep='')) ) {
      # reduce the number of cells programmed to prolferate
      assign( paste('gen',k,'p',sep='') , get(paste('gen',k,'p',sep='')) - 1 )
      # determine the fate of the 2 daughter cells
      U = sum(runif(2, min=0, max=1) < P_rand[4])
      if (U == 2) { # both daughter cell are programmed to proliferate
        assign(paste('gen',k+1,'p',sep = '') , get(paste('gen',k+1,'p',sep = '')) + 2 )
      } 
      else if(U == 0) { # both daughter cells are  programmed to die
        assign(paste('gen',k+1,'d',sep = '') , get(paste('gen',k+1,'d',sep = '')) + 2 )
      }
      else{ # one is programmed to proliferate and one is programmed to die with prob 2P[death]P[div]
        assign(paste('gen',k+1,'p',sep = '') , get(paste('gen',k+1,'p',sep = '')) + 1 )
        assign(paste('gen',k+1,'d',sep = '') ,  get(paste('gen',k+1,'d',sep = '')) + 1 )
      }
      
    }
    
    if(delta_t %in% get(paste('li.d',k,sep=''))) {
      assign(paste('gen',k,'d',sep = '') ,  get(paste('gen',k,'d',sep = '')) - 1 )
      
    } 
  }else{
    
    if ( (x%% 2) == 0 ) { # check if the minimum index is even or odd to determine position
      k = (0.5*x) - 1
    } else{               # do this for odd position 
      k = (0.5*x) - 0.5   
    }
    # cells that go onto the last generation poliferate into division classes of no interest in this study 
    if ( delta_t %in% li.p20){
      gen20p = gen20p - 1  # reduce the number of cells programmed to prolferate
    } 
    if (delta_t %in% li.d20){
      gen20d = gen20d - 1  # reduce the number of cells programmed to die
    }
    
  }
  # we append cells within 0.001 of the time points.
  if (is.element(signif(Tt, 1), t.pts)) { # measurements at time points saved in csv file
    if (abs(Tt - tpts) < 0.001 ) {
      N_prol <- c(gen0p ,gen1p , gen2p, gen3p ,gen4p , gen5p, gen6p ,gen7p , gen8p,
                  gen9p ,gen10p , gen11p, gen12p ,gen13p , gen14p, gen15p ,gen16p ,
                  gen17p,gen18p ,gen19p , gen20p )
      N_death <- c(gen0d ,gen1d , gen2d, gen3d ,gen4d , gen5d, gen6d ,gen7d , gen8d,
                   gen9d ,gen10d , gen11d, gen12d ,gen13d , gen14d, gen15d ,gen16d ,
                   gen17d ,gen18d ,gen19d , gen20d )
      
     
      Ncell[1,1] <- N_init
      
      for (t in q:q) {  
        for (n in 1:21 ) {
          Ncell[t,n] <- N_prol[n] + N_death[n]
        }
      }
      q <- q + 1
      tpts <- tpts + 12 
    }
      
   
  }
  
  Tt = Tt + delta_t       # update the time. loop breaks when sum of time is > Tmax
  #print(Tt)               # uncomment to visualize time points
  
}







