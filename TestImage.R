#library("jpeg")
library("magick")


# Download image
y = image_read("https://cdn.arstechnica.net/wp-content/uploads/2016/02/5718897981_10faa45ac3_b-640x624.jpg")
image_info(y)

#
plot(y)
image_trim(y,0)
image_scale(y, "150x150")

plot(y)




