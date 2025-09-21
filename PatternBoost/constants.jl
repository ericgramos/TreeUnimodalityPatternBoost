const NB_LOCAL_SEARCHES::Int64 = 120 # number of local searches to do at once. 
                                      # probably best to set it divisible by number of threads, if using multithreading
                                      # increase this number if episodes are too fast

const NUM_INITIAL_EMPTY_OBJECTS = 50_000 #on the first run, this many rollouts will be done from empty graphs/objects
const FINAL_DATABASE_SIZE = 50_000 # Size of learning set
const TARGET_DB_SIZE = 5_000 # Database size will never be twice bigger than that
const OBJ_TYPE = Vector{Int}  # type for vectors encoding constructions
const REWARD_TYPE = Float32  # type for rewards 
