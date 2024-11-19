function run_ICG(input_file, output_file)
load(input_file);
res = ICG(X);
save(output_file,'res');
end
