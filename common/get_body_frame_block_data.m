function [rin_body, rout_body] = get_body_frame_block_data(rin_block, rout_block, q_block, R_block)
%GET_BODY_FRAME_BLOCK_DATA Map one particle's source and target nodes to body coordinates.

rin_body = (R_block' * (rin_block - q_block)')';
rout_body = (R_block' * (rout_block - q_block)')';

end
