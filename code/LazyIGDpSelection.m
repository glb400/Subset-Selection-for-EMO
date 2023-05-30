 function [Subset,res,time] = LazyIGDpSelection(PopObj,selNum)
    tic;
    a = size(PopObj, 1);
    M = size(PopObj, 2);
    selectedPop = [];
    IGDArray = zeros(1, a);

    heap.size=0;
    heap.index=[];
    heap.hvc=[];
    
    for n = 1 : a
       heap.index(n) = n; 
       currentPop = [selectedPop; PopObj(n, :)];
       [~, heap.hvc(n)] = IGDCp(currentPop, PopObj, size(currentPop, 1), IGDArray);
       heap.size = heap.size +1;  
        if toc > 3600 % The max time for one run is one hour
            break
        end   
    end
    toc

    si = heap.size;
    while(si >= 1)
        heap = heap_down_sink(si, heap);
        si = si - 1;
    end
    %% Select first solution
    top = heap.index(1);
    selectedPop = [selectedPop; PopObj(top, :)];
    [IGDArray, ~] = IGDCp(selectedPop, PopObj, size(selectedPop, 1), IGDArray);
    heap = min_heap_popup(heap);
    count = 0;
    
    %% Select other solutions
    while size(selectedPop, 1) < selNum
        while true
            last_index = heap.index(1);
            top = heap.index(1);
            cuurentPop = [selectedPop; PopObj(top, :)];
            [~, heap.hvc(1)] = IGDCp(cuurentPop, PopObj, size(cuurentPop, 1), IGDArray);
            heap = heap_down_sink(1, heap);
            if(heap.index(1) == last_index)            
                selectedPop = [selectedPop; PopObj(heap.index(1), :)];
                [IGDArray, ~] = IGDCp(selectedPop, PopObj, size(selectedPop, 1), IGDArray);
                heap = min_heap_popup(heap);
      
                break;  
            end
            if toc > 3600 % The max time for one run is one hour
                break
            end 
        end
        if toc > 3600 % The max time for one run is one hour
            break
        end 
    end
    
    %% output 
    Subset = selectedPop;
    time = toc;
    res = IGDp(Subset,PopObj);
    disp(res);
end

function heap = heap_down_sink(self, heap)
lchild =self*2;
rchild =self*2+1;
if(lchild <= heap.size && rchild <= heap.size)
    if(heap.hvc(self) < heap.hvc(lchild) || heap.hvc(self) < heap.hvc(rchild))
        if(heap.hvc(lchild) >= heap.hvc(rchild))
            
            heap = heap_swap_node(lchild,self,heap);
            heap = heap_down_sink(lchild, heap);
        else
            heap = heap_swap_node(rchild,self, heap);
            heap = heap_down_sink(rchild, heap);
        end
    end
elseif(lchild <= heap.size)
    if(heap.hvc(self) < heap.hvc(lchild))
        heap = heap_swap_node(lchild,self, heap);
    end
end
end

function heap = heap_swap_node(node1,node2, heap)

    t_index = heap.index(node1);
    t_hvc = heap.hvc(node1);
    
    heap.index(node1) = heap.index(node2);
    heap.hvc(node1) = heap.hvc(node2);

    heap.index(node2) = t_index;
    heap.hvc(node2) = t_hvc;

end


function heap = min_heap_popup(heap)
heap = heap_swap_node(1,heap.size, heap);
heap.size =heap.size -1;
heap = heap_down_sink(1,heap);
end

