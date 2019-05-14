### Setup

SignalSetCommand <- setClass(
  Class = 'SignalSetCommand',
  slots = c(
    name = 'character',
    time.stamp = 'POSIXct',
    call.string = 'character',
    params = 'ANY'
  )
)

NetSet <- setClass(Class ="NetSet",

                   slots = c(
                     graph_data = 'data.frame',
                     graph_object = 'list'
                   )
)

GroupSet <- setClass(Class = "GroupSet",

                     slots = c(
                       group = "character",
                       pD = 'data.frame',
                       fD = 'data.frame',
                       DE = "list",
                       AggBulk = "matrix",
                       network = "NetSet",
                       commands = "SignalSetCommand"),

                     prototype = list(group = "ALL",
                                      pD = as.data.frame(matrix(NA, nrow = 50, ncol =2)),
                                      fD = as.data.frame(matrix(NA, nrow = 100, ncol =2)),
                                      DE = vector(mode = "list"),
                                      AggBulk = matrix(NA, nrow = 100, ncol = 2),
                                      network = new('NetSet'),
                                      commands = new('SignalSetCommand'))
)

SignalSet  <- setClass(Class = "SignalSet",
                       slots = c(
                         counts = "dgCMatrix",
                         norm_counts = "dgCMatrix",
                         group_ids = "character", # A vector of the possible groups
                         group_membership = "matrix", # A vector, length of all cells, that tells us does it belong to a group or not
                         group_sets = "list", # A list of group sets, on for each group_id
                         active_group = "GroupSet"
                       ),

                       prototype = list(counts = Matrix::Matrix(c(0,1), nrow = 100, ncol = 50, sparse = T),
                                        norm_counts = Matrix::Matrix(c(0,1), nrow = 100, ncol = 50, sparse = T),
                                        group_ids = c("ALL"),
                                        group_membership = matrix(c(rep(1,50))),
                                        group_sets = list(new("GroupSet")),
                                        active_group = new("GroupSet"))

)

# Need to write a creator!! Takes in a matrix!


# First we initialize the SignalSet, and this creates a group, ALL
# We then do some intial processing
# From this we id potential "groups"
# Need a create_group function
# Need a active_group function
# Need a write_group function? When we are working it is saving into active group. IF we then set another active group, what happens to changes?
# Need to decide,


