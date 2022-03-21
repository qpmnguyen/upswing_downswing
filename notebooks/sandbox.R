
h <- settings$h
x <- df_inc %>% pull(time_value)
y <- df_inc %>% pull(cases)
x_ref <- x[6]
ii <- x > x_ref - h & x <= x_ref + h
print(x_ref)
print(x[ii])
xx <- as.numeric(x[ii])
x_ref <- as.numeric(x_ref)
yy <- y[ii]
right = xx > x_ref
left = xx <= x_ref
b = mean(yy[right])
a = mean(yy[left])
hh = mean(xx[right]) - mean(xx[left])
(b/a - 1)



surge_cum %>% filter(time_value == x[6] + 4) %>% pull(rel_change)
surge_inc %>% filter(time_value == x[6]) %>% pull(growth_adj)
