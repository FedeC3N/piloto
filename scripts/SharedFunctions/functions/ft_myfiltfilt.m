function data = ft_myfiltfilt ( b, a, data )

data = ft_checkdata ( data, 'datatype', 'raw' );

for trial = 1: numel ( data.trial )
    data.trial { trial } = myfiltfilt ( b, a, data.trial { trial }' )';
end
