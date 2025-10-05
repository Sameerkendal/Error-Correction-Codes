
corrupted_mp3_path = 'music8_corrupted.mp3';
output_mp3_path ='Group8_decoded.mp3';
%% Parameters
n = 125;
k = 100;
tc = floor((n - k)/2);         
prim_poly = 285;              

%% Read the corrupted file (raw bytes)
fid = fopen(corrupted_mp3_path, 'rb');
if fid < 0
    error('Cannot open input file: %s', corrupted_mp3_path); 
end

bytes = fread(fid, inf, 'uint8=>uint8');
fclose(fid);


orig_size = 64000;
if numel(bytes) < orig_size
    warning('Input file has only %d bytes (< %d). Will decode what is present.', numel(bytes), orig_size);
end
bytes = bytes(:);  


num_blocks = ceil(numel(bytes)/n);
pad_needed = num_blocks*n - numel(bytes);
if pad_needed > 0
    bytes = [bytes; zeros(pad_needed,1,'uint8')];
end


evalpts = gf(0:n-1, 8, prim_poly);    

%% Output buffer 

out_bytes = zeros(num_blocks*k, 1, 'uint8');

%% Process each RS block of n bytes

write_idx = 1;
disp(num_blocks);
for b = 1:num_blocks
    disp(b);
    block = bytes((b-1)*n + (1:n));         
    y = gf(block.', 8, prim_poly);  

    decoded_msg = berlekamp_welch_decoder(evalpts,y,n, k, tc, prim_poly);
    disp(size(decoded_msg,2));
    
    out_bytes(write_idx : write_idx + k - 1) = uint8(decoded_msg.x);
    write_idx = write_idx + k;
end


out_bytes = out_bytes(1:orig_size);

%% Write the decoded MP3
fidw = fopen(output_mp3_path, 'wb');
if fidw < 0
    error('Cannot open output file: %s', output_mp3_path); 
end
fwrite(fidw, out_bytes, 'uint8');
fclose(fidw);

fprintf('Decoding complete. Wrote %d bytes to %s\n', numel(out_bytes), output_mp3_path);






