function invf(fY)
    Y = fY.^3
    i = (Y .< 0.008856)
    Y[i] = (fY[i] .- 4/29)*(108/841)
    return Y
end

function gammacorrection(R)
    Rp = real(1.099*R.^0.45 .- 0.099)
    i = (R .< 0.018)
    Rp[i] = 4.5138 * R[i]
    return Rp
end


function invgammacorrection(Rp)
    R = real(((Rp .+ 0.099)./1.099).^(1/0.45))
    i = (R .< 0.018)
    R[i] = Rp[i] / 4.5138;
    return R;
end


function f(Y)
    fY = real((Y.+0im).^(1/3));
    i = (Y .< 0.008856);
    fY[i] = Y[i] * (841/108) .+ (4/29);
    return fY;
end

function hsv(Image,SrcSpace)

    Image = rgb(Image,SrcSpace);
    V = max(Image,[],3);
    S = (V .- min(Image,[],3))./(V .+ (V .== 0));

    Image[:,:,1] = rgbtohue(Image);
    Image[:,:,2] = S;
    Image[:,:,3] = V;
    return Image;
end

function huetorgb(m0, m2, H)
    N = size(H);
    H = min.(max.(H[:],0),360)/60;
    m0 = m0[:];
    m2 = m2[:];
    F = H .- round.(H ./ 2) .* 2;
    M = [m0, m0 .+ (m2 .- m0) .* abs.(F), m2];
    Num = length(m0);
    j = [2 1 0;1 2 0;0 2 1;0 1 2;1 0 2;2 0 1;2 1 0] * Num;
    k = Int64.(floor.(H)) .+ 1;

    col1 = Int64.(((j[k,1]' .+ reshape.([1:Num]',1,Num)[1])[1] - 1) / 4 + 1)
    col2 = Int64.(((j[k,2]' .+ reshape.([1:Num]',1,Num)[1])[1] - 1) / 4 + 1)
    col3 = Int64.(((j[k,3]' .+ reshape.([1:Num]',1,Num)[1])[1] - 1) / 4 + 1)
    Image = [reshape(M[col3], Num, 1), reshape(M[col2], Num, 1), reshape(M[col1], Num, 1)]';
    return Image;
end

function hsl(Image,SrcSpace)
    if SrcSpace == "hsv"
        MaxVal = Image[:,:,3];
        MinVal = (1 .- Image[:,:,2]).* MaxVal;
        L = 0.5*(MaxVal .+ MinVal);
        temp = min.(L,1 .- L);
        Image[:,:,2] = 0.5*(MaxVal .- MinVal) ./ (temp .+ (temp .== 0));
        Image[:,:,3] = L;
    else
       Image = rgb(Image,SrcSpace);
       MinVal = minimum(Image, dims=3);
       MaxVal = maximum(Image, dims=3);
       L = 0.5*(MaxVal .+ MinVal);
       temp = min.(L,1 .- L);
       S = 0.5*(MaxVal .- MinVal)./(temp .+ (temp == 0));
       Image[:,:,1] = rgbtohue(Image);
       Image[:,:,2] = S;
       Image[:,:,3] = L;
    end

    return Image;
end

function lab(Image,SrcSpace)
    WhitePoint = [0.950456,1,1.088754];

    if SrcSpace == "lab"
       return Image;
    elseif SrcSpace == "lch"
       C = Image[:,:,2];
       Image[:,:,2] = cos.(Image[:,:,3]*pi./180).*C;
       Image[:,:,3] = sin.(Image[:,:,3]*pi./180).*C;
    else
       Image = xyz(Image,SrcSpace);

       X = Image[:,:,1]./WhitePoint[1];
       Y = Image[:,:,2]./WhitePoint[2];
       Z = Image[:,:,3]./WhitePoint[3];
       fX = f(X);
       fY = f(Y);
       fZ = f(Z);
       Image[:,:,1] = 116*fY .- 16;
       Image[:,:,2] = 500*(fX .- fY);
       Image[:,:,3] = 200*(fY .- fZ);
    end

    return Image;
end

function luv(Image,SrcSpace)

    WhitePoint = [0.950456,1,1.088754];
    WhitePointU = (4*WhitePoint[1])./(WhitePoint[1] + 15*WhitePoint[2] + 3*WhitePoint[3]);
    WhitePointV = (9*WhitePoint[2])./(WhitePoint[1] + 15*WhitePoint[2] + 3*WhitePoint[3]);

    Image = xyz(Image,SrcSpace);
    U = (4*Image[:,:,1])./(Image[:,:,1] + 15*Image[:,:,2] + 3*Image[:,:,3]);
    V = (9*Image[:,:,2])./(Image[:,:,1] + 15*Image[:,:,2] + 3*Image[:,:,3]);
    Y = Image[:,:,2]/WhitePoint[2];
    L = 116*f(Y) .- 16;
    Image[:,:,1] = L;
    Image[:,:,2] = 13*L.*(U .- WhitePointU);
    Image[:,:,3] = 13*L.*(V .- WhitePointV);
    return Image;
end


function lch(Image,SrcSpace)
    Image = lab(Image,SrcSpace);
    H = atan.(Image[:,:,3], Image[:,:,2]);
    H = H*180/pi + 360*(H .< 0);
    Image[:,:,2] = sqrt.(Image[:,:,2].^2 + Image[:,:,3].^2);
    Image[:,:,3] = H;
    return Image;
end

function rgbtohue(Image)
    M = mapslices(sort, Image, dims = 3)
    i = mapslices(sortperm, Image, dims = 3)
    i = i[:,:,3]
    Delta = M[:, :, 3] - M[:, :, 1]
    Delta = Delta + (Delta .== 0);
    R = Image[:, :, 1]
    G = Image[:, :, 2]
    B = Image[:, :, 3]
    H = zeros(size(R))
    k = (i .== 1)
    k_index = findall(x->x==true, k[:])
    H[k_index] = (G[k_index] - B[k_index]) ./ Delta[k_index]
    k = (i .== 2)
    k_index = findall(x->x==true, k[:])
    H[k_index] = 2 .+ (B[k_index] - R[k_index]) ./ Delta[k_index]
    k = (i .== 3)
    k_index = findall(x->x==true, k[:])
    H[k_index] = 4 .+ (R[k_index] - G[k_index]) ./ Delta[k_index]
    H = 60 .* H + 360 .* (H .< 0)
    Delta_index = findall(x->x==0, Delta[:])
    H[Delta_index] .= NaN
    return H
end


function rgb(Image,SrcSpace)
    if SrcSpace == "rgb"
       return Image;
    elseif SrcSpace == "hsv"
       Image = huetorgb((1 .- Image[:,:,2]) .* Image[:,:,3],Image[:,:,3],Image[:,:,1]);
    elseif SrcSpace == "hsl"
       L = Image[:,:,2];
       Delta = Image[:,:,2].* min.(L,1 .- L);
       Image = huetorgb(L-Delta,L+Delta,Image[:,:,1]);
    elseif SrcSpace in ["xyz", "lab", "luv", "lch"]

       Image = xyz(Image,SrcSpace);

       T = [3.240479 -1.53715 -0.498535; -0.969256 1.875992 0.041556; 0.055648 -0.204043 1.057311]
       R = T[1]*Image[:,:,1] + T[4]*Image[:,:,2] + T[7]*Image[:,:,3];
       G = T[2]*Image[:,:,1] + T[5]*Image[:,:,2] + T[8]*Image[:,:,3];
       B = T[3]*Image[:,:,1]+ T[6]*Image[:,:,2] + T[9]*Image[:,:,3];

       AddWhite = -min.(min.(min.(R,G),B),0);
       Scale = max.(max.(max.(R,G),B) .+ AddWhite,1);
       R = (R .+ AddWhite) ./ Scale;
       G = (G .+ AddWhite) ./ Scale;
       B = (B .+ AddWhite) ./ Scale;

       Image[:,:,1] = gammacorrection(R)
       Image[:,:,2] = gammacorrection(G)
       Image[:,:,3] = gammacorrection(B)
    else
       T = gettransform(SrcSpace);
       temp = inv(T[:,1:3]);
       T = [temp,-temp.*T[:,4]];
       R = T[1]*Image[:,:,1] + T[4]*Image[:,:,2] + T[7]*Image[:,:,3] + T[10]
       G = T[2]*Image[:,:,1] + T[5]*Image[:,:,2] + T[8]*Image[:,:,3] + T[11]
       B = T[3]*Image[:,:,1]+ T[6]*Image[:,:,2] + T[9]*Image[:,:,3] + T[12]
       AddWhite = -min.(min.(min.(R,G),B),0);
       Scale = max.(max.(max.(R,G),B)+AddWhite,1);
       R = (R .+ AddWhite)./Scale;
       G = (G .+ AddWhite)./Scale;
       B = (B .+ AddWhite)./Scale;
       Image[:,:,1] = R;
       Image[:,:,2] = G;
       Image[:,:,3] = B;
    end

    Image = min.(max.(Image,0),1);
    return Image;
end


function xyz(Image,SrcSpace)
    # Convert to CIE XYZ from 'SrcSpace'
    WhitePoint = [0.950456 1 1.088754]
    if SrcSpace == "xyz"
        return Image
    elseif SrcSpace == "luv"
        #Convert CIE L*uv to XYZ
        WhitePointU = (4 * WhitePoint[1]) ./ (WhitePoint[1] + 15 * WhitePoint[2] + 3 * WhitePoint[3])
        WhitePointV = (9 * WhitePoint[2]) ./ (WhitePoint[1] + 15 * WhitePoint[2] + 3 * WhitePoint[3])
        L = Image[:, :, 1]
        Y = (L .+ 16) ./ 116
        Y = invf(Y) * WhitePoint[2];
        U = Image[:, :, 2] ./ (13 .* L + 1e-6 .* (L .== 0)) .+ WhitePointU
        V = Image[:, :, 3] ./ (13 .* L + 1e-6 .* (L .== 0)) .+ WhitePointV
        Image[:, :, 1] = -(9 .* Y .* U) ./ ((U .- 4) .* V - U .* V) #X
        Image[:, :, 2] = Y #Y
        Image[:, :, 3] = (9 .* Y - (15 .* V .* Y) - (V .* Image[:, :, 1])) ./ (3 .* V)  #Z
    elseif (SrcSpace == "lab" || SrcSpace == "lch")
        Image = lab(Image,SrcSpace)
        # Convert CIE L*ab to XYZ
        fY = (Image[:, :, 1] .+ 16) ./ 116
        fX = fY + Image[:,:,2] ./ 500
        fZ = fY - Image[:,:,3] ./ 200
        Image[:, :, 1] = WhitePoint[1] .* invf(fX) #X
        Image[:, :, 2] = WhitePoint[2] .* invf(fY) #Y
        Image[:, :, 3] = WhitePoint[3] .* invf(fZ)
    else
        #Convert from some gamma-corrected space
        #Convert to Rec. 701 R'G'B'
        Image = rgb(Image,SrcSpace)
        #Undo gamma correction
        R = invgammacorrection(Image[:, :, 1])
        G = invgammacorrection(Image[:, :, 2])
        B = invgammacorrection(Image[:, :, 3])
        #Convert RGB to XYZ
        T = inv([3.240479 -1.53715 -0.498535; -0.969256 1.875992 0.041556; 0.055648 -0.204043 1.057311])
        Image = float.(Image)
        Image[:,:,1] = T[1] .* R + T[4] .* G + T[7] .* B  # X
        Image[:,:,2] = T[2] .* R + T[5] .* G + T[8] .* B  # Y
        Image[:,:,3] = T[3] .* R + T[6] .* G + T[9] .* B  # Z
    end
    return Image
end

function gettransform(Space)
    # Get a colorspace transform: either a matrix describing an affine transform,
    # or a string referring to a conversion subroutine
    if Space == "ypbpr"
        T = [0.299 0.587 0.114 0; -0.1687367 -0.331264 0.5 0; 0.5 -0.418688 -0.081312 0]
    elseif Space == "yuv"
        # R'G'B' to NTSC/PAL YUV
        # Wikipedia: http://en.wikipedia.org/wiki/YUV
        T = [0.299 0.587 0.114 0; -0.147 -0.289 0.436 0; 0.615 -0.515 -0.100 0]
    elseif Space == "ydbdr"
        #R'G'B' to SECAM YDbDr
        #Wikipedia: http://en.wikipedia.org/wiki/YDbDr
        T = [0.299 0.587 0.114 0; -0.450 -0.883 1.333 0; -1.333 1.116 0.217 0]
    elseif Space == "yiq"
        #R'G'B' in [0,1] to NTSC YIQ in [0,1];[-0.595716,0.595716];[-0.522591,0.522591]
        #Wikipedia: http://en.wikipedia.org/wiki/YIQ
        T = [0.299 0.587 0.114 0; 0.595716 -0.274453 -0.321263 0; 0.211456 -0.522591 0.311135 0]
    elseif Space == "ycbcr"
        # R'G'B' (range [0,1]) to ITU-R BRT.601 (CCIR 601) Y'CbCr
        # Wikipedia: http://en.wikipedia.org/wiki/YCbCr
        # Poynton, Equation 3, scaling of R'G'B to Y'PbPr conversion
        T = [65.481 128.553 24.966 16; -37.797 -74.203 112.0 128; 112.0 -93.786 -18.214 128]
    elseif Space == "jpegycbcr"
        #Wikipedia: http://en.wikipedia.org/wiki/YCbCr
        T = [0.299 0.587 0.114 0; -0.168736 -0.331264 0.5 0.5; 0.5 -0.418688 -0.081312 0.5] .* 255
    elseif (Space == "rgb" || Space == "xyz" || Space == "hsv" || Space == "hsl" || Space == "lab" || Space == "luv" || Space == "lch")
        T = Space
    else
        error("Unknown color space, ", Space, ".")
    end
    return T
end

function alias(Space)
    Space = replace(Space, "cie"=>"")
    if isempty(Space)
        Space = "rgb"
    end
    if (Space == "ycbcr" || Space == "ycc")
        Space = "ycbcr"
    elseif (Space == "hsv" || Space == "hsb")
        Space = "hsv"
    elseif (Space == "hsl" || Space == "hsi" || Space == "hls")
        Space = "hsl"
    elseif (Space == "rgb" || Space == "yuv" || Space == "yiq" || Space == "ydbdr" || Space == "ycbcr" || Space == "jpegycbcr" || Space == "xyz" || Space == "lab" || Space == "luv" || Space == "lch")
        return Space
    end
    return Space
end

function Parse(Str)
    Conversion = Str
    # Parse conversion argument
    if typeof(Str) == String
        Str = lowercase(replace(replace(Str,"-"=>"")," "=>""))
        k = findall(">", Str)
        if length(k) == 1         #Interpret the form 'src->dest'
            k = k[1][1]
            SrcSpace = Str[1:k-1]
            DestSpace = Str[k+1:end]
        else
            k = findall("<", Str)
            if length(k) == 1      #Interpret the form 'dest<-src'
                k = k[1][1]
                DestSpace = Str[1:k-1]
                SrcSpace = Str[k+1:end];
            else
                error("Invalid conversion, ", Str, ".")
            end
        end
        SrcSpace = alias(SrcSpace)
        DestSpace = alias(DestSpace)
    else
        SrcSpace = 1       # No source pre-transform
        DestSpace = Conversion
        if any(size(Conversion) .!= 3)
            error("Transformation matrix must be 3x3.")
        end
    end
    return SrcSpace, DestSpace
end

function colorspace(Conversion,varargin)

    SrcSpace, DestSpace = Parse(Conversion)
    Image = varargin

    FlipDims = (size(Image,3) == 1)


    if FlipDims
        Image = permutedims(Image,[1,3,2])
    end

    if eltype(Image) != Float64
        Image = Float64.(Image) ./ 255
    end

    if size(Image,3) != 3
        error("Invalid input size.")
    end


    SrcT = gettransform(SrcSpace);
    DestT = gettransform(DestSpace);

    if ~all(isletter, SrcT) & ~all(isletter, DestT)
        T = [DestT[:,1:3]*SrcT[:,1:3],DestT[:,1:3]*SrcT[:,4]+DestT[:,4]]
        Temp = zeros(size(Image))
        Temp[:,:,1] = T[1].*Image[:,:,1] .+ T[4]*Image[:,:,2] .+ T[7]*Image[:,:,3] .+ T[10]
        Temp[:,:,2] = T[2].*Image[:,:,1] .+ T[5]*Image[:,:,2] .+ T[8]*Image[:,:,3] .+ T[11]
        Temp[:,:,3] = T[3].*Image[:,:,1] .+ T[6]*Image[:,:,2] .+ T[9]*Image[:,:,3] .+ T[12]
        Image = Temp;
    elseif ~all(isletter, DestT)
        Image = rgb(Image,SrcSpace)
        Temp = zeros(size(Image))
        Temp[:,:,1] = DestT[1].*Image[:,:,1] .+ DestT[4].*Image[:,:,2] .+ DestT[7].*Image[:,:,3] .+ DestT[10]
        Temp[:,:,2] = DestT[2].*Image[:,:,1] .+ DestT[5].*Image[:,:,2] .+ DestT[8].*Image[:,:,3] .+ DestT[11]
        Temp[:,:,3] = DestT[3].*Image[:,:,1] .+ DestT[6].*Image[:,:,2] .+ DestT[9].*Image[:,:,3] .+ DestT[12]
        Image = Temp;
    else
        Image = eval(Symbol(DestT))(Image, SrcSpace)
    end

    if FlipDims
        Image = permute(Image,[1,3,2])
    end
    varargout = Image;

    return varargout
end
