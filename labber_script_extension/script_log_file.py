import numpy as np
import sys
import os
import copy
sys.path.append(os.environ.get(
    'LABBERPATH', r'C:\Program Files (x86)\Labber\Script'))  # NOQA: E402
import Labber as Labber


class ScriptLogFile(Labber.LogFile):
    """
        Extension of Labber LogFile that provides some additional features for reading the data.
    """
    def yieldData(self, log_channel=None, inner=0):
        """
        Returns the data at the index.as a numpy array.

        Args:
            log_channel (str): name of the log channel which data is returned.
            inner (int): number of inner loop dimensions to yield in chunks.
        Returns:
            The data as a float or numpy array.
        """
        (channel_values, channel_names, dimensions) = self.getDataMatrixMetaData(log_channel)

        if inner >= len(dimensions):
            yield self.getDataMatrix(log_channel)[0]
        else:

            n_elements = np.prod(dimensions[:inner])

            if channel_names[0] == 'vector_data' and False:
                loop_start = 0
                if inner == 0:
                    n_elements = 1
                else:
                    n_elements = np.prod(dimensions[:inner])
            else:
                loop_start = 1
                n_elements = np.prod(dimensions[loop_start:inner+1])
                if inner > 0:
                    n_elements = n_elements//dimensions[inner]

            for ind in range(np.prod(dimensions[max(inner, loop_start):])):
                data = []
                for ii in range(n_elements):
                    temp_data = self.getData(log_channel, entry=ind*n_elements + ii)
                    if inner == 0:
                        for jj in range(dimensions[0]):
                            yield temp_data[jj]
                    else:
                        data.extend(temp_data)
                data = np.array(data)
                if len(dimensions[:inner] > 0):
                    dimensions_reversed = dimensions[:inner]
                    dimensions_reversed = dimensions_reversed[::-1]
                    data = np.reshape(data, dimensions_reversed)
                    data = np.transpose(data)
                if inner != 0:
                    yield data

    def getDataMatrixMetaData(self, log_channel=None):
        """
        Returns the metadata for getDataMatrix.

        Args:
            log_channel (str): name of the log channel which information is returned.
        Returns:
            A tuple with three elements. The first contains a list of arrays with the step channel data, the second element
            contains a list of step channel names, and the third element contains a tuple of the dimensions of the data.
        """

        step_channels = self.getStepChannels()
        log_channels = self.getLogChannels()
        if log_channel is None:
            log_channel = log_channels[0]['name']

        log_channel_dict = None
        for log_channel_dict_iter in log_channels:
            if log_channel_dict_iter['name'] == log_channel:
                log_channel_dict = log_channel_dict_iter
                break
        channel_names = []
        channel_values = []
        dimensions = []
        for channel in step_channels:
            if isinstance(channel['values'], np.ndarray) and len(channel['values']) > 1:
                if 'unit' in channel and channel['unit'] != '':
                    channel_names.append(channel['name'] + ' (' + channel['unit'] + ')')
                else:
                    channel_names.append(channel['name'])
                channel_values.append(channel['values'])
                dimensions.append(len(channel['values']))

        # XXX log is not specified. Check if this is intentional
        data = self.getData(log_channel, entry=0)
        if log_channel_dict['vector']:
            dimensions.insert(0, data.shape[-1])
            channel_names.insert(0, 'vector_data')
            channel_values.insert(0, np.linspace(0, 1, data.shape[-1]))

        # XXX in getDataMatrix there is check if len(dimensions) == 0. Should probably be here, too.
        first_element = dimensions.pop(0)
        dimensions = dimensions[::-1]
        dimensions.append(first_element)
        dimensions = np.array(dimensions[::-1])

        if 'unit' in log_channel_dict and log_channel_dict['unit'] != '':
            channel_names.append(log_channel + ' (' + log_channel_dict['unit'] + ')')
        else:
            channel_names.append(log_channel)

        return (channel_values, channel_names, dimensions)

    def getDataMatrix(self, log_channel=None, log=-1, channel_order=None, take_index=None):
        """
        Returns the data as a numpy array.

        Args:
            log_channel (str): name of the log channel which data is returned. If not specified, the first one is returned.
            log (int): Index of log which data to retrieve. Defaults to the last log (-1)
            channel_order (list): list of step_channels in order. The data is reordered so that the first dimension corresponds to the first step channel in the list.
                                  If step channel corresponds to a singleton dimension, a new dimension is created in its place. As a result.
            take_index (int or tuple of ints): Only used if step_channels is not None. If there are more non-singleton step channels than specified by step_channels, the other dimensions are picked by take_index. As a result, the number of returned dimensions is equal to number of specified step_channels.
        Returns:
            A tuple with three elements. The first is an ndarray with the data, second contains a list of arrays with the step channel data, and the third element
            contains a list of step channel names and the log channel name as the last element.
        """
        step_channels = self.getStepChannels()
        log_channels = self.getLogChannels()
        if log_channel is None:
            log_channel = log_channels[0]['name']

        log_channel_dict = None
        for log_channel_dict_iter in log_channels:
            if log_channel_dict_iter['name'] == log_channel:
                log_channel_dict = log_channel_dict_iter
                break
        channel_names = []
        channel_names_wo_units = []
        channel_values = []
        dimensions = []
        for channel in step_channels:
            if isinstance(channel['values'], np.ndarray) and len(channel['values']) > 1:
                if 'unit' in channel and channel['unit'] != '':
                    channel_names.append(channel['name'] + ' (' + channel['unit'] + ')')
                else:
                    channel_names.append(channel['name'])
                channel_names_wo_units.append(channel['name'])
                channel_values.append(channel['values'])
                dimensions.append(len(channel['values']))

        data = self.getData(log_channel, log=log)
        if log_channel_dict['vector']:
            dimensions.insert(0, data.shape[-1])
            channel_names.insert(0, 'vector_data')
            channel_names_wo_units.insert(0, 'vector_data')
            channel_values.insert(0, np.linspace(0, 1, data.shape[-1]))

        if len(dimensions) == 0:
            return (np.squeeze(data), channel_values, channel_names)
        first_element = dimensions.pop(0)
        dimensions = dimensions[::-1]
        dimensions.append(first_element)

        # if data.size != np.prod(dimensions):
        #     missing_factor = np.prod(dimensions) / data.size
        #     data = data[:int(dimensions[0]*missing_factor), :]
        #     dimensions
        data = np.reshape(data, dimensions)
        data = np.transpose(data)
        if 'unit' in log_channel_dict and log_channel_dict['unit'] != '':
            channel_names.append(log_channel + ' (' + log_channel_dict['unit'] + ')')
        else:
            channel_names.append(log_channel)

        #print(channel_names)
        #print(channel_names_wo_units)

        current_order = []
        new_channels = []
        #new_channels_with_units = []
        new_channel_values = []
        channel_dict = self.getChannelValuesAsDict()
        if channel_order is not None:
            for channel in channel_order:
                try:
                    current_order.append(channel_names_wo_units.index(channel))
                except ValueError:
                    # channel is not in the list
                    if channel in channel_dict:
                        new_channels.append(channel)
                        new_channel_values.append([channel_dict[channel]])
                        current_order.append(len(channel_names_wo_units) + len(new_channels) - 1)
                        # XXX Deal with adding the units
                    else:
                        print(f'channel {channel:s} is not in the list.')
                        raise
                #print(current_order[-1])

        new_order = list(np.arange(0, len(channel_names_wo_units) + len(new_channels)))

        def diff(first, second):
            second = set(second)
            return [item for item in first if item not in second]

        new_order = current_order + diff(new_order, current_order)

        if len(new_channels) > 0:
            data = np.expand_dims(data, axis=list(np.arange(-len(new_channels), 0)))
        data = np.transpose(data, axes=new_order)
        channel_names_wo_units.extend(new_channels)
        channel_names_wo_units = [channel_names_wo_units[ii] for ii in new_order]
        channel_values.extend(new_channel_values)
        channel_values = [channel_values[ii] for ii in new_order]
        channel_names = channel_names_wo_units

        if take_index is not None and channel_order is not None:
            if len(channel_values) > len(channel_order):
                data = np.transpose(data)
                data = data[take_index]
                data = np.transpose(data)
                channel_values = channel_values[:len(channel_order)]
                channel_names = channel_names[:len(channel_order)]

        channel_names.append(log_channel)

        return (data, channel_values, channel_names)
