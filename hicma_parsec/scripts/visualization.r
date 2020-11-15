library('ggplot2')



starsh = read.csv('starsh_gen.csv')
starsh = starsh[starsh$type == 8, ] # they approximate_L task, the key 8 is obtained from h5 event_types table
starsh$time = starsh$end - starsh$begin
p = ggplot(starsh, aes(n, m)) + geom_tile(aes(fill = time), colour = "white") + 
                                scale_fill_gradient(low = "white", high = "steelblue")
# output the plot
#p
