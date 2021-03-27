# input:
# image: input image
# r: radius of window ([x-r:x+r, y-r:y+r])

# Output:
# Output(x, y)=sum(sum(imSrc(x-r:x+r,y-r:y+r)))
function boxfilter(image, r)
    h, w = size(image)
    out_image = zeros(size(image))
    for i=1:h
        for j=1:w
            x1 = max(min(i-r, h), 1)
            x2 = max(min(i+r, h), 1)
            y1 = max(min(j-r, w), 1)
            y2 = max(min(j+r, w), 1)
            out_image[i, j] = sum(image[x1:x2, y1:y2])
        end
    end
    return out_image
end
