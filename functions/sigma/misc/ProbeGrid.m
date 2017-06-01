classdef ProbeGrid < handle
    properties(SetAccess=private)
        % main properties
        iProbe
        jProbe
        maxSize
        % derived properties
        fields
        numValues
    end
    
    methods
        function obj = ProbeGrid(iProbe, jProbe, maxSize, fieldNames)
        %%
        % Class constructor.
        % 
        % Written by P. DERIAN 2017-02-05.
            obj.maxSize = maxSize;
            obj.iProbe = iProbe;
            obj.jProbe = jProbe;
            init = NaN([numel(iProbe), numel(jProbe), maxSize]);
            for i=1:numel(fieldNames)
                obj.fields.(fieldNames{i}) = init;
                obj.numValues.(fieldNames{i}) = 0;
            end
        end
        
        function update_field(self, fieldName, data)
        %%
        % Update given probe field with fresh data.
        %
        % Written by P. DERIAN 2014-02-05.
            kProbe = self.numValues.(fieldName) + 1;
            if kProbe>self.maxSize
                warning('ProbeGrid:update_field:maxSizeReached',...
                        'maximum probe size (%d) of field %s reached, this and subsequent update_field() ignored.',...
                        self.maxSize, fieldName);
            end
            self.fields.(fieldName)(:,:,kProbe) = data(self.iProbe, self.jProbe);
            self.numValues.(fieldName) = kProbe;
        end
        
        function data = get_field(self, fieldName)
        %%
        % Retrieve a field.
        %
        % Written by P. DERIAN 2017-02-05.
            data = self.fields.(fieldName)(:,:,1:self.numValues.(fieldName));
        end
        
        function obj = saveobj(self)
            obj.iProbe = self.iProbe;
            obj.jProbe = self.jProbe;
            obj.maxSize = self.maxSize;
            obj.fields = self.fields;
            obj.numValues = self.numValues;
        end
    end
    
    methods(Static)
        function obj = loadobj(s)
            if isstruct(s)
                newobj = ProbeGrid(s.iProbe, s.jProbe, s.maxSize, ...
                                   fieldnames(s.fields));
                newobj.fields = s.fields;
                newobj.numValues = s.numValues;
                obj = newobj;
            else
                obj = s;
            end
        end
    end
end