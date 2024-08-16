function matrix_output = Matrix_shuffle(matrix_input)

    % for a 2d matrix, randomly changes its order.
    % version 1
    %     [x,y] = size(matrix_input);
    %     matrix_temp = reshape(matrix_input,x*y,1);
    %     order_temp = randperm(x*y);
    %     matrix_temp = matrix_temp(order_temp);
    %     matrix_output = reshape(matrix_temp,[x,y]);
        
    % version 2
    % [x,y] = size(matrix_input);
    % order_temp = reshape(randperm(x*y,x*y),x,y);
    % matrix_output = matrix_input(order_temp);
    
    % version 3
    [x,y] = size(matrix_input);
    order_temp = reshape(Shuffle2(1:(x*y)),x,y);
    matrix_output = matrix_input(order_temp);
end