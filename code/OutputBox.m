classdef OutputBox < handle
    %Esta clase representará el historial de información que se le irá
    %proveiendo al usuario.
    
    properties
        text;
    end
    
    methods
        %Constructor
        function obj = OutputBox(texto)
            obj.text = texto;
        end
        
        %Función para añadir texto 
        function text = addText(obj,texto)
            if(texto ~= " ")
                obj.text = sprintf('%s\n\n%s', "> " + texto) + obj.text;
                text = obj.text;
            else
                obj.text = obj.text + sprintf('%s\n\n%s',texto);
                text = obj.text;
            end
             
        end
       
    end
    
end

