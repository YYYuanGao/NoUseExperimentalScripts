%��ֹʵ�����
function abort  
    ShowCursor;              %��ʾ���
    Screen('CloseAll');      %�ر�������
    %reset_test_gamma;       %ȥ��gamma����
    warning on;              %���������ʾwarning
    error('aborting due to user cancellation...');
end