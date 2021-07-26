from collections import Iterable
import datetime
import logging
import numpy as np
import sys
import os
import copy
sys.path.insert(0, os.environ.get(
    'LABBERPATH', r'C:\Program Files\Labber\Script'))   # NOQA: E402
import Labber as Labber
from Labber import ScriptTools


class ScriptObject(ScriptTools.MeasurementObject, Labber.Scenario):
    """
    A Labber measurement object based on Labber.Scenario. Contains also features from ScriptTools.MeasurementObject.
    """

    def __init__(self, sCfgFileIn=None, sCfgFileOut=None):
        Labber.Scenario.__init__(self, sCfgFileIn)
        #self.file_in = sCfgFileIn
        self.file_out = sCfgFileOut
        #self.log_name = sCfgFileOut
        #self.scenario = ScriptTools.load_scenario_as_dict(self.file_in)
        #self.scenario = Labber.Scenario(self.file_in)
        self.primary_channel = None  # XXX

    # For backwards compability. Scenario object uses log_name
    @property
    def file_out(self):
        return self.log_name

    @file_out.setter
    def file_out(self, s):
        self.log_name = s

    def performMeasurement(self, return_data=True, **kwargs):
        """ Performs the measurement.

        Args:
            return_data (bool): boolean indicating whether data should be returned by the labber measurement.
            **kwargs (dict): Possible keywords: use_scheduler (bool). This is supported only by labber versions.
        """
        temp_file_name = os.path.splitext(self.file_out)[0]+'_tmp.labber'
        file_dir_path = os.path.dirname(self.file_out)

        # Create subdirectories if they don't exist.
        if file_dir_path.strip():
            os.makedirs(file_dir_path, exist_ok=True)

        #ScriptTools.save_scenario_as_binary(self.scenario, temp_file_name)
        self.save(temp_file_name)
        labber_meas_object = ScriptTools.MeasurementObject(
            temp_file_name, self.file_out)
        if self.primary_channel is not None:
            if hasattr(labber_meas_object, 'setPrimaryChannel'):
                labber_meas_object.setPrimaryChannel(self.primary_channel)
            else:
                labber_meas_object.setMasterChannel(self.primary_channel)
        #return labber_meas_object.performMeasurement(return_data, use_scheduler)  # use_scheduler is not supported by the older versions of labber
        return labber_meas_object.performMeasurement(return_data, **kwargs)

    def save_as_binary(self, filename):
        """ Saves the scipt object as labber binary (.labber file).

        Args:
            filename (str): name of the saved binary.
        """
        #ScriptTools.save_scenario_as_binary(self.scenario, filename)
        self.save(filename)

    def setMasterChannel(self, channel_name):
        """
            Sets the master channel.

            When other parameters are updated, creates a lookup table for them with respect to the master channel.

            Does not work correctly at the moment. This is depracated.
        """

        logging.warning('setMasterChannel() is depracated. Calling setPrimaryChannel instead.')
        self.setPrimaryChannel(channel_name)
        # raise Exception('setMasterChannel not implementd.')
        #self.master_channel = channel_name
        # temp_file_name = os.path.splitext(self.file_out)[0]+'_tmp.labber'
        # ScriptTools.save_scenario_as_binary(self.scenario, temp_file_name)
        # labber_meas_object = ScriptTools.MeasurementObject(temp_file_name,self.file_out)
        # labber_meas_object.setMasterChannel(channel_name)

    def setPrimaryChannel(self, channel_name):
        """
            Sets the master channel.

            When other parameters are updated, creates a lookup table for them with respect to the master channel.

        Args:
            channel_name (str): Name of the primary channel
        """
        self.primary_channel = channel_name

    def setOutputFile(self, filename):
        """
        Sets the name of the output file.
        """
        self.file_out = filename

    def rearrangeLog(self, channel_name, *extra_arg):
        """
        This is not implemented, yet.
        """
        raise Exception('Not implemented.')

    def updateInstrumentValue(self, instrument_name, parameter_name, value, full_name=None):
        """
        Updates a value directly in the instrument (corresponds to changing values in the left-hand side in labber)

        Arguments
        instrument_name (str) -- name of the instrument
        parameter_name (str) -- name of the parameter
        value (float) -- value of the new parameter
        full_name (str) -- Full labber name of the channel ('[instrument_name] - [parameter_name]'). If given, overrides instrument_name and parametr_name
        """

        if full_name is not None:
            split_name = full_name.split(' - ', 1)
            instrument_name = split_name[0]
            parameter_name = split_name[1]

        for current_instrument in self.instruments:
            if instrument_name == current_instrument.com_config.name:
                current_instrument.values[parameter_name] = value
                return

        logging.warning('Instrument ' + instrument_name + ' not found.')
        return

    def copy_instrument_values(self, source_instrument_name, target_instrument_name):
        """ Copies the instrument values, except for communication related parameters. Instruments have to be of same type.
        """
        source = self.get_instrument(source_instrument_name)
        target = self.get_instrument(target_instrument_name)
        if source is not None and target is not None:
            target.values = copy.deepcopy(source.values)

    # Replaced by Scenario's function
    # def get_instrument(self, instrument_name):
    #     """ Returns the instrument with given name.

    #     """
    #     for current_instrument in self.scenario['instruments']:
    #         if instrument_name == current_instrument['com_config']['name']:
    #             return current_instrument
    #     return None

    def removeStepChannel(self, channel):
        """ Removes the step channel if it exists.

        Args:
            channel (str): channel to be removed
        """
        if channel in self.step_names():
            self.remove_step(channel)

    def getStepChannels(self):
        """
        Returns all the step channels (right side of the measurement program).
        """
        return self.step_items

    def getStepChannelNames(self, filter_string=None, instrument_name=None):
        """
        Returns all the step channels which name matches the filter string.

        Arguments:
            filter_string (str) -- string which is used to filter the returned step channels.
            instrument_name (str) -- if not None, only stepChannels that are related to the given instrument string are returned.

        Returns:
            A list of step channel names that match the provided filter strings.
        """
        step_channels = self.getStepChannels()
        step_channel_names = []

        for channel in step_channels:
            channel_name = channel.channel_name

            filter_passed = False
            instrument_name_passed = False

            if filter_string is None or filter_string in channel_name:
                filter_passed = True

            split_name = channel_name.split(' - ', 1)
            if len(split_name) > 1:
                instrument_name_part = split_name[0]
                parameter_name = split_name[1]

                if instrument_name is None or instrument_name_part == instrument_name:
                    instrument_name_passed = True
            else:
                parameter_name = split_name[0]
                if instrument_name is None:
                    instrument_name_passed = True

            if filter_passed and instrument_name_passed:
                step_channel_names.append(channel_name)

        return step_channel_names

    def getChannels(self):
        """
        Returns all the channels (left side of the measurement program).

        """
        return self.channels

    def getChannel(self, name):
        """
        Get a reference to channel by its name.

        Arguments:
        name (str) -- Name of the channel.

        Returns:
        Reference to the channel or None if not found.
        """

        # if name in self.getChannelNames():
        # for channel in self.getChannels():
        #     if channel['instrument'] + ' - ' + channel['quantity'] == name:
        #         return channel
        #     elif 'name' in channel.keys() and channel['name'] == name:
        #         return channel
        # return None
        if name in self.channel_names():
            return self.get_channel(name)
        else:
            return None

    def getStepChannel(self, name):
        """ Get a reference to a step channel by its name.

        Args:
            name (str): Name of the step channel.

        Returns:
            Reference to the step channel or None if not found.
        """
        for ii, step_channel_name in enumerate(self.getStepChannelNames()):
            if step_channel_name == name:
                return self.getStepChannels()[ii]

    def getChannelNames(self):
        """
        Returns a list of all channel names.
        """
        # channel_names = []
        # for channel in self.getChannels():
        #     if 'name' in channel.keys():
        #         channel_names.append(channel['name'])
        #     else:
        #         channel_names.append(
        #             channel['instrument'] + ' - ' + channel['quantity'])
        # return channel_names
        return self.channel_names()

    def getInstrumentValueNames(self, filter_string=None, instrument_name=None, is_case_sensitive=False):
        """
        Returns the names of the instrument values that match the filter strings.

        Arguments:
            filter_string (str) -- String which is used to filter the returned instrument value names.
            instrument_name (str) -- If not None, only instrument value names that are related to the given instrument string are returned.
            is_case_sensitive (bool, default=False) -- Boolean indicating whether filter matching is case sensitive.

        Returns:
            A list of instrument value names that match the provided filter strings.
        """

        if not is_case_sensitive and filter_string is not None:
            filter_string = filter_string.lower()
        if not is_case_sensitive and instrument_name is not None:
            instrument_name = instrument_name.lower()

        instruments = self.instruments
        # for current_instrument in self.scenario['instruments']:
        instrument_value_names = []
        for instrument in instruments:
            current_instrument_name = instrument.com_config.name
            for value_name in instrument.values.keys():
                filter_passed = False
                instrument_name_passed = False

                if not is_case_sensitive:
                    value_name_to_test = value_name.lower()
                if filter_string is None or filter_string in value_name_to_test:
                    filter_passed = True
                if not is_case_sensitive:
                    instrument_name_to_test = current_instrument_name.lower()

                if instrument_name is None or instrument_name == instrument_name_to_test:
                    instrument_name_passed = True
                if filter_passed and instrument_name_passed:
                    instrument_value_names.append(
                        current_instrument_name + ' - ' + value_name)

        return instrument_value_names

    def getLogChannels(self):
        """
        Returns a list of all log channels in the scenario.
        """
        #log_channels = self.scenario['log_channels']
        log_channels = self.log_channels
        return log_channels

    def getRangeItemValues(self, range_item):
        """ Returns values represented by range_item.

        Args:
            range_item (Labber.config.step.RangeItem): Labber range object
        Returns:
            Numpy array of values represented by range item.
        """
        if range_item.range_type == Labber.config.step.RangeType.SINGLE:
            return [range_item.single]
        elif range_item.range_type == Labber.config.step.RangeType.CENTERSPAN:
            if range_item.step_type == Labber.config.step.RangeStep.N_PTS:
                return range_item.center + np.linspace(-range_item.span/2, range_item.span/2, range_item.n_pts)
            else:
                return range_item.center + np.arange(-range_item.span/2, range_item.span/2+range_item.step, range_item.step)
        elif range_item.range_type == Labber.config.step.RangeType.STARTSTOP:
            if range_item.step_type == Labber.config.step.RangeStep.N_PTS:
                return np.linspace(range_item.start, range_item.stop, range_item.n_pts)
            else:
                return np.arange(range_item.start, range_item.stop+range_item.step, range_item.step)

        logging.warning('range_item has invalued value of RangeType or RangeStep')
        return None

    def getStepChannelValues(self, channel_name):
        """ Returns the values of a step channel in a list.

        Note that this function does not currently support lookup tables or equations, and returns an incorrect value if they are used.

        Args:
            channel_name (str): Name of the step_channel.

        Returns:
            Numpy array of channel values.
        """
        channel = self.getStepChannel(channel_name)
        if channel is None:
            logging.warning(f'Step channel {channel_name} not found. Returning None.')
            return None
        values = np.array([])
        for range_item in channel.range_items:
            values = np.append(values, self.getRangeItemValues(range_item))
        return values

    def getChannelValue(self, channel_name):
        """ Returns the value of a channel.

        Args:
            channel_name (str): Name of the channel which value is returned

        Returns:
            value of the channel or None is channel is not found.

        """
        split_name = channel_name.split(' - ', 1)
        if len(split_name) > 1:
            instrument_name = split_name[0]
            quantity_name = split_name[1]
        else:
            # Channel Name is a nick name, since no dash was found.
            channel = self.getChannel(channel_name)
            if channel is not None:
                instrument_name = channel.instrument
                quantity_name = channel.quantity
            else:
                return None
        for current_instrument in self.instruments:
            if instrument_name == current_instrument.com_config.name:
                if quantity_name in current_instrument.values:
                    return current_instrument.values[quantity_name]
                else:
                    return None

    def printStepChannels(self, filter_string=None, instrument_name=None, verbose=False):
        """
        Prints all the step channels which name matches the filter string.

        Arguments:
            filter_string(str) - - string which is used to filter the returned step channels.
            instrument_name(str) - - if not None, only stepChannels that are related to the given instrument string are returned.
            verbose(bool) - - if True, all available information about the channel is listed.
        """
        step_channels = self.getStepChannels()

        step_channel_names = []
        step_channel_data = []

        for channel in step_channels:
            channel_name = channel.channel_name

            filter_passed = False
            instrument_name_passed = False

            if filter_string is None or filter_string in channel_name:
                filter_passed = True

            split_name = channel_name.split(' - ', 1)
            if len(split_name) > 1:
                instrument_name_part = split_name[0]
                parameter_name = split_name[1]

                if instrument_name is None or instrument_name_part == instrument_name:
                    instrument_name_passed = True
            else:
                parameter_name = split_name[0]
                if instrument_name is None:
                    instrument_name_passed = True

            if filter_passed and instrument_name_passed:
                step_channel_names.append(channel_name)
                step_channel_data.append(channel)

        if not verbose:
            print(*step_channel_names, sep="\n")
        else:
            for ii in range(len(step_channel_names)):
                pass

    def move_channel_to(self, channel_name, index):
        """
        Moves a step channel to the given index in the list of step channels.

        This affects the order in witch the paramters are set when initializing the instruments.
        Can also be used the change the order in which step channels are looped. The first channel is the innermost loop.

        Arguments:
            channel_name(str) - - Name of the step channel
            index(int) - - new position of the channel in the list of channels. If negative, added to the end of stepchannels. Note that -1 adds to the end (different than most functions in python).
        """

        step_channels = self.getStepChannels()

        if channel_name in self.getStepChannelNames():
            if index < 0:
                index = len(step_channels) + index  # Note the additional +1. -1 moves to the end.
            self.set_step_position(channel_name, index)
        # for ii in range(len(step_channels)):
        #     channel = step_channels[ii]
        #     if channel_name == channel.channel_name:
        #         step_channels.insert(index, channel)
        #         if(ii < index):
        #             step_channels.pop(ii)
        #         else:
        #             step_channels.pop(ii+1)
        #         return

    def move_log_channel_to(self, channel_name, index):
        """
        Moves a log channel to the given index in the list of log channels.

        Arguments:
            channel_name(str) - - Name of the log channel
            index(int) - - new position of the channel in the list of channels.
        """

        log_channels = self.getLogChannels()

        if channel_name in log_channels:
            if index < 0:
                index = len(log_channels) + index  # Note the additional +1. -1 moves to the end.
            self.set_log_position(channel_name, index)

    def add_equation(self, channel_name, equation, variables):
        """ A helper function to add channel relation.

        Args:
            channel_name (str): Name of the channel for which the equation is added.
            equation (str): The equation string. Ex: 'p1 + 2'
            variables (dict of strs): A dictionary specifying the variables in the equation. Ex: {'p1': 'Instrument Name - Channel Name', 'p2': 'Instrument Name 2 - Another Channel'}. You can create lookup table in the following way:
        {'p1': {channel_name: 'Instrument Name - Channel Name', lookup_x: [list of values], lookup_y: [list of values]}, 'p2': 'Instrument Name 2 - Another Channel'}
        """
        self.updateStepChannelsByDict({channel_name: {'EQ': equation, 'VARS': variables}})

        # Make sure that the step channels exist
        for key, value in variables.items():
            if type(value) is dict:
                # Lookup table, add required channels. NOT IMPLEMENTED.
                pass
            elif self.getStepChannel(value) is None and value != 'Step values':
                self.addStepChannel(value)

    def write_lookup(self, lookup_table, channel_name, ref_channel_name, add_current_value=True, interp='Linear'):
        """ A wrapper to write a lookup table for a channel in equation. Not that this will override the existing equation. The full functionality can be achieved with add_equation() function.

    Args:
        lookup_table (list): Lookup table as list in Labber convention.
        channel_name (Str): Channel name of channel that obtains lookup table.
        ref_channel_name (Str): Reference channel name.
        add_current_value (bool, optional): If True, current value of the channel is added to the lookup table.
        interp (str): Interpolation function. Can be one of ['Linear', 'Nearest', 'Quadratic', 'Cubic'].
    """

        if add_current_value:
            self.updateStepChannelsByDict({
                channel_name: {
                    'EQ': 'x+p1',
                    'VARS': {
                        'x': 'Step values', 'p1': {'channel_name': ref_channel_name, 'lookup_x': lookup_table[0], 'lookup_y': lookup_table[1], 'interp': interp}
                    }
                }
            })
        else:
            self.updateStepChannelsByDict({channel_name: {'EQ': 'p1', 'VARS': {'x': 'Step values',
                'p1': {'channel_name': ref_channel_name, 'lookup_x': lookup_table[0], 'lookup_y': lookup_table[1], 'interp': interp}}}})

    def update_step_parameters(self, channel_name,
                               sweep_rate=None,
                               wait_after=None,
                               after_last=None,
                               final_value=None,
                               sweep_mode=None,
                               use_outside_sweep_rate=None,
                               sweep_rate_outside=None,
                               alternate_direction=None):
        """ A wrapper that can be used to update channel step parameters.

        Args:
            channel_name (str): Name of the step channel. If it doesn't exist, the channel is created.
            sweep_rate (double): Rate at which the value should be changed (unit/s)
            wait_after (double): Number of seconds to wait after the value is set.
            after_last (str): What to do after sweep is finished. Either 'Goto first point', 'Stay at final', or 'Goto value...' (yes, there are three dots in the string).
            sweep_mode (str): Must be one of ['Between points', 'Off', 'Continuous'].
            final_value (double): Final value for the instrument. If given, will set after_last to 'Goto value...'.
            use_outside_sweep_rate (bool): Boolean indicating whether to use different sweep rate outside.
            sweep_rate_outside (double): Sweep rate outside of loop. If value other than zero is given, use_outside_sweep_rate is set to True.
            alternate_direction (bool): Boolean indicating whether to sweep in alternate_direction
        """
        param_dict = {}
        if wait_after is not None:
            param_dict['wait_after'] = wait_after

        if final_value is not None:
            param_dict['final_value'] = final_value
            if final_value:
                param_dict['after_last'] = 'Goto value...'
        if after_last is not None:
            if after_last in ['Goto first point', 'Stay at final', 'Goto value...']:
                param_dict['after_last'] = after_last
            else:
                logging.warning('after_last must be on of "Goto first point", "Stay at final", "Goto value..."')
        if sweep_mode is not None:
            if sweep_mode in ['Between points', 'Off', 'Continuous']:
                param_dict['sweep_mode'] = sweep_mode
            else:
                logging.warning('Sweep mode must be on of "Between points", "Off", or "Continuous"')

        if sweep_rate_outside is not None:
            if sweep_rate_outside:
                param_dict['use_outside_sweep_rate'] = True
            param_dict['sweep_rate_outside'] = sweep_rate_outside
        if use_outside_sweep_rate is not None:
            param_dict['use_outside_sweep_rate'] = use_outside_sweep_rate
        if alternate_direction is not None:
            param_dict['alternate_direction'] = alternate_direction

        self.updateStepChannelsByDict({channel_name: {'PARAM': param_dict}})

        if sweep_rate is not None:
            channel = self.getStepChannel(channel_name)
            if channel is not None and channel['sweep_mode'].lower() != 'off':
                for step_item in channel['step_items']:
                    step_item['sweep_rate'] = sweep_rate

    def updateValue(self, channel_name, value, itemType='SINGLE', step_index=None):
        """
        Update a value in the config file.

        Arguments:
            channel_name(str) - - Name of channel to update.
            value(float) - - New value to set  to channel.
            itemType(str, optional) - - Step item parameter to set, must be one
                                        of {SINGLE, START, STOP, CENTER, SPAN,
                                        STEP, N_PTS, PARAM, VARS, EQ}. Default is SINGLE.
                                        If PARAM type is chosen, value is a dict
                                        of name-value pairs from
                                         'wait_after': 0.0,
                                         'final_value': 0.0,
                                         'show_advanced': False,
                                         'use_relations': False,
                                         'equation': 'x',
                                         'step_unit': 'Instrument',
                                         'after_last': 'Goto first point',
                                         'sweep_mode': 'Off',
                                         'use_outside_sweep_rate': False,
                                         'sweep_rate_outside': 0.0,
                                         'alternate_direction': False

                                         EQ is used to set the equation for channel relations. Value should be a string representing the equation.
                                         VARS is used to define the variables used in the the channel relation equations. The value is a dict, which keys are the
                                         variable names used in the equations and the values are the channel names to which the variables refer.
                                         Alternatively, the variables can be defined to be lookup tables.

                                         In order to use channel relations, the option has to be turned on using PARAM use_relations: True.

            step_index(int, optional) - - index in the step items which value will be updated. If None, all other steps will be cleared.
        """
        if value is None:
            logging.warning('Trying to update channel ' + channel_name + ' with None value')

        def set_sweep_parameter(channel, value, itemType, index):
            # if(not np.isscalar(value)):
            #    channel['relation_parameters'][0]['use_lookup'] = True
            #    channel['relation_parameters'][0]['lookup'] = value
            #    return
            if(itemType == 'PARAM'):
                for param_name, param_value in value.items():
                    setattr(channel, param_name, param_value)
                    #channel[param_name] = param_value
                return
            if(itemType == 'EQ'):
                #channel['equation'] = value
                setattr(channel, 'equation', value)
                set_sweep_parameter(channel, {
                    'show_advanced': True,
                    'use_relations': True
                }, 'PARAM', index)
                return
            if(itemType == 'VARS'):
                # Constructs relation_parameters object
                #setattr(channel, 'relation_parameters', [])
                relation_parameter_list = []
                if 'x' not in value.keys():
                    relation_parameter = Labber.config.step.RelationParameter()
                    relation_parameter.variable = 'x'
                    relation_parameter.channel_name = 'Step values'
                    relation_parameter.lookup = None
                    relation_parameter_list.append(relation_parameter)
                for dict_param, dict_value in value.items():
                    relation_parameter = Labber.config.step.RelationParameter()
                    relation_parameter.variable = dict_param
                    relation_parameter.lookup = Labber.config.step.LookUpTable()
                    if isinstance(dict_value, dict):
                        relation_parameter.use_lookup = True
                        # if 'interp' in dict_value.keys():
                        #     interp = dict_value['interp']
                        # else:
                        #     interp = 'Linear'
                        for relation_key, relation_value in dict_value.items():
                            if relation_key == 'lookup_x':
                                relation_parameter.lookup.xdata = np.array(dict_value['lookup_x'])
                            elif relation_key == 'lookup_y':
                                relation_parameter.lookup.ydata = np.array(dict_value['lookup_y'])
                            elif relation_key == 'interp':
                                relation_parameter.lookup.interp = relation_value
                            else:
                                setattr(relation_parameter, relation_key, relation_value)
                    else:
                        relation_parameter.channel_name = dict_value
                        relation_parameter.lookup = None
                    relation_parameter_list.append(relation_parameter)
                set_sweep_parameter(channel, {
                    'show_advanced': True,
                    'use_relations': True}, 'PARAM', index)
                channel.relation_parameters = relation_parameter_list
                return

            if index is None:
                if (len(channel.range_items)) > 0:
                    del channel.range_items[1:]
                index = 0

            if(index >= len(channel.range_items)):
                for ii in range(index-len(channel.range_items)+1):
                    range_item = Labber.config.step.RangeItem()
                    range_item.n_pts = 1
                    channel.range_items.append(range_item)
                    # channel.range_items.append({
                    #     'range_type': 'Start - Stop',
                    #     'step_type': 'Fixed # of pts',
                    #     'single': 0,
                    #     'start': 0,
                    #     'stop': 1,
                    #     'center': 0.5,
                    #     'span': 1,
                    #     'step': 0.1,
                    #     'n_pts': 1,
                    #     'interp': 'Linear',
                    #     'sweep_rate': 0.0
                    # })
            if(itemType == 'SINGLE'):
                channel.range_items[index].range_type = 'Single'
                channel.range_items[index].single = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'START'):
                channel.range_items[index].range_type = 'Start - Stop'
                channel.range_items[index].start = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'STOP'):
                channel.range_items[index].range_type = 'Start - Stop'
                channel.range_items[index].stop = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)

            if(itemType == 'SPAN'):
                channel.range_items[index].range_type = 'Center - Span'
                channel.range_items[index].span = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'CENTER'):
                channel.range_items[index].range_type = 'Center - Span'
                channel.range_items[index].center = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'N_PTS'):
                channel.range_items[index].step_type = 'Fixed # of pts'
                channel.range_items[index].n_pts = int(value)
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'STEP'):
                channel.range_items[index].step_type = 'Fixed step'
                channel.range_items[index].step = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)

        step_channels = self.getStepChannels()
        for current_channel in step_channels:
            # print(current_channel['channel_name'],channel_name)
            if(channel_name == current_channel.channel_name):
                # und match: ' + channel_name + ' and ' + current_channel['channel_name'])
                set_sweep_parameter(current_channel, value,
                                    itemType, step_index)
                # Extend here with advanced options!
                return
        # channel not found, add a new channel.
        new_channel = self.add_step(channel_name)
        #new_channel = Labber.config.step.StepItem()
        # new_channel = {'channel_name': channel_name,
        #                'wait_after': 0.0,
        #                'final_value': 0.0,
        #                'show_advanced': False,
        #                'use_relations': False,
        #                'equation': 'x',
        #                'step_unit': 'Instrument',
        #                'after_last': 'Goto first point',
        #                'sweep_mode': 'Off',
        #                'use_outside_sweep_rate': False,
        #                'sweep_rate_outside': 0.0,
        #                'alternate_direction': False,
        #                'step_items': [{
        #                    'range_type': 'Single',
        #                    'step_type': 'Fixed # of pts',
        #                    'single': 0,
        #                    'start': 0,
        #                    'stop': 1,
        #                    'center': 0.5,
        #                    'span': 1,
        #                    'step': 0.1,
        #                    'n_pts': 11,
        #                    'interp': 'Linear',
        #                    'sweep_rate': 0.0
        #                }],
        #                'relation_parameters': [{
        #                    'variable': 'x',
        #                    'channel_name': 'Step values',
        #                    'use_lookup': False, 'lookup': None}],
        #                'optimizer_config': {
        #                    'Enabled': False,
        #                    'Start value': 7660000000.0,
        #                    'Initial step size': 1532000000.0,
        #                    'Min value': 7660000000.0,
        #                    'Max value': 15320000000.0,
        #                    'Precision': 766000.0}
        #                }
        #step_channels.append(new_channel)
        # print(channel_name)

        # In order to designate the channel to the correct instrument, instrument name needs to be extacted.
        # The following approach fails, if the is ' - ' in the instrument name or there is no instrument name at all.
        # This can be the case if nick name is used for the channel. However, in that case it should always exist as a step channel.
        if channel_name is None:
            logging.warning('Trying to setup step channel, but the channel name is None.')
        labber_channel = self.getChannel(channel_name)
        if labber_channel is None:
            # Channel is not found. Need to create one.
            split_name = channel_name.split(' - ', 1)

            if len(split_name) > 1:
                instrument_name = split_name[0]
                parameter_name = split_name[1]
                self.addChannel(instrument_name, parameter_name)
            else:
                # Channel Name is a nick name, since no dash was found.
                raise ValueError('Parameter ' + channel_name +
                                 ' your tried to update does not exist and cannot be created. Try creating it first using addChannel().')

        set_sweep_parameter(new_channel, value, itemType, step_index)

        logging.info('Step channel ' + channel_name + ' not found. It has been created.')
        return

    def setSignalConnection(self, target, source):
        """
        Connects two signals.

        Only signals with the same types can be used as target - source pairs. This is not checked. Currently does not work if the channel has nick name.

        Arguments:
        target - Name of the channel that uses the source signal.
        source - Name of the channel that acts as a source signal for the target.
        """

        target_channel = self.getChannel(target)
        source_channel = self.getChannel(source)

        if target_channel is None:
            self.addChannel(None, None, full_name=target)
            target_channel = self.getChannel(target)
        if source_channel is None:
            self.addChannel(None, None, full_name=source)
            source_channel = self.getChannel(source)

        self.add_connection(source, target)
        #target_channel.signal_source = source

    '''
    def removeSignalConnection(self, target):
        """Deletes a signal connection.

        Args:
            target (str): Target channel for which the signal connection is to be removed.
        """
        del(self.getChannel(target)['signal_source'])
    '''

    def setSignalConnectionsByDict(self, connection_dict):
        """
        Creates signal connections from dictionary of target - source pairs.

        Arguments:
            connection_dict(dict) - - keys are the target channel names, values are the source channel names.
        """
        for target, source in connection_dict.items():
            self.setSignalConnection(target, source)

    def addStepChannel(self, step_channel_name, value=None):
        """
        Adds a new step channel with default parameters.

        Does not check if channel exists in channels. If a step channel of same name already exists, new channel is not added.

        Parameters:
            step_channel_name(str) - - The name of the step channel.
            value (double), default(None): the initial value of the step channel. If None, the value of the original channel is used (if it exists).
        """
        if value is None:
            value = self.getChannelValue(step_channel_name)
        if value is None:
            value = 0.0
        # new_channel = {'channel_name': step_channel_name,
        #                'wait_after': 0.0,
        #                'final_value': 0.0,
        #                'show_advanced': False,
        #                'use_relations': False,
        #                'equation': 'x',
        #                'step_unit': 'Instrument',
        #                'after_last': 'Goto first point',
        #                'sweep_mode': 'Off',
        #                'use_outside_sweep_rate': False,
        #                'sweep_rate_outside': 0.0,
        #                'alternate_direction': False,
        #                'step_items': [{
        #                    'range_type': 'Single',
        #                    'step_type': 'Fixed # of pts',
        #                    'single': value,
        #                    'start': 0,
        #                    'stop': 1,
        #                    'center': 0.5,
        #                    'span': 1,
        #                    'step': 0.1,
        #                    'n_pts': 11,
        #                    'interp': 'Linear',
        #                    'sweep_rate': 0.0
        #                }],
        #                'relation_parameters': [{
        #                    'variable': 'x',
        #                    'channel_name': 'Step values',
        #                    'use_lookup': False, 'lookup': None}],
        #                'optimizer_config': {
        #                    'Enabled': False,
        #                    'Start value': 7660000000.0,
        #                    'Initial step size': 1532000000.0,
        #                    'Min value': 7660000000.0,
        #                    'Max value': 15320000000.0,
        #                    'Precision': 766000.0}
        #                }

        if step_channel_name in self.getStepChannelNames():
            pass
        else:
            self.add_step(step_channel_name)

        if value is not None:
            self.updateValue(step_channel_name, value)

        # step_channels = self.getStepChannels()
        # step_channels.append(new_channel)

    def copy_step_channel(self, old_channel_name, new_channel_name):
        old_step_channel = self.getStepChannel(old_channel_name)
        if old_step_channel is not None:
            new_step_channel = copy.deepcopy(old_step_channel)
            new_step_channel.channel_name = new_channel_name
            self.removeStepChannel(new_channel_name)
            self.getStepChannels().append(new_step_channel)

    def addChannel(self, instrument_name, quantity, name=None, full_name=None, physical_unit=None, instrument_unit=None, value=None):
        """
        Adds a new channel to the measurement object.

        Arguments:
            instrument_name(str): name of the instrument the channel belongs to.
            quantity(str): the name of the channel in the instrument driver
            name(str), default(None): the name which the channel is called
            full_name(str), default(None): The full name of the channel. If not None, this overrides insturment_name and quantity.
            physical_unit (str), default(None): The physical unit of the channel (Hz etc). If None, no unit is added.
            instrument_unit (str), default(None): The insturment unit of the channel (Hz etc). If None, no unit is added.
            value (double), default(None): Updates the value of the channel. If None, the default value of the instrument is used, or if channel already exists,
                                           its value is left unchanged..
        """

        if full_name is not None:
            split_name = full_name.split(' - ', 1)
            instrument_name = split_name[0]
            quantity = split_name[1]

        if value is not None:
            self.updateInstrumentValue(instrument_name, quantity, value)

        for channel in self.channels:
            if(channel.instrument == instrument_name and channel.quantity == quantity):
                # Channel already exitst
                if name is not None:
                    self.nameChannel(instrument_name + ' - ' + quantity, name)
                if physical_unit is not None:
                    channel.unit_physical = physical_unit
                if instrument_unit is not None:
                    channel.unit_instrument = instrument_unit
                return

        channel = Labber.config.scenario.Channel()
        channel.instrument = instrument_name
        channel.quantity = quantity
        # channel = {
        #     'instrument': instrument_name,
        #     'quantity': quantity,
        # }
        if physical_unit is not None:
            channel.unit_physical = physical_unit
        if instrument_unit is not None:
            channel.unit_instrument = instrument_unit
        if name is not None:
            channel.name = name
        self.channels.append(channel)

    def nameChannel(self, old_name, new_name):
        """ Gives a nick name to the channel.

        In addition to naming the channel, renames all the step channels and log channels to which
        old_name of the channel refers to. Currently does not rename signal connections.
        Args:
            old_name (str): Full name of the channel in format 'instrument - quantity' or the nick name of the channel.
            nick_name (str): The name with which the channel should be called. If None, nick Name is removed.
        """
        channel = self.getChannel(old_name)

        if channel is not None:
            full_channel_name = channel.instrument + ' - ' + channel.quantity
            if channel.name is not None:
                nick_name = channel.name
            else:
                nick_name = None
            if new_name is not None:
                channel.name = new_name
            else:
                if hasattr(channel, 'name'):
                    channel.name = None
            step_channel = self.getStepChannel(old_name)
            if step_channel is None:
                step_channel = self.getStepChannel(nick_name)
            if step_channel is not None:
                if new_name is not None:
                    step_channel.channel_name = new_name
                else:
                    step_channel.channel_name = full_channel_name
            for ii, log_channel in enumerate(self.getLogChannels()):
                if log_channel == old_name or log_channel == nick_name:
                    if new_name is not None:
                        self.getLogChannels()[ii] = new_name
                    else:
                        self.getLogChannels()[ii] = channel.instrument + ' - ' + channel.quantity

    def setParameter(self, name, value):
        """
        Sets a parameter value pair for the scenario.

        Currently available name/value pairs are
        'log_parallel': bool - - ??
        'log_name': str - - name of the log file to be created
        'comment': str - - comment of the log
        'wait_between': float - - ??
        'time_per_point': float - - should not be touched
        'arm_trig_mode': bool - - ??
        'hardware_loop': bool - - enable hardware looping
        'trig_channel': str - - name of the trigger channel in the form "Instrument - channel name"
        'logger_mode': bool - - ??
        """

        # TODO somehow test the validity of the given parameters, e.g check if the instruments support hardware looping
        logging.warning('This function is depracated in Scenario object. Most of these items are now in settings.')
        if hasattr(self, name):
            setattr(self, name, value)
        else:
            logging.warning('Parameter ' + name + ' does not exist.')

    def setComment(self, comment):
        """
        Add a comment to the measurement.
        """
        if comment is not None:
            self.comment = comment

    def setTags(self, tags, override=True):
        """
        Sets the tags for the measurement object.

        Arguments:
                tags(list of str or str) - - List tags to be set or a single tag
                override(bool) - - If True, all the existing tags will be overridden. If False, appends tags.
        """
        if isinstance(tags, list):

            if override:
                self.tags.tags = tags
            else:
                self.tags.tags.extend(tags)
        else:
            if override:
                self.tags.tags = [tags]
            else:
                self.tags.tags.append(tags)

    def setUser(self, user):
        """
        Sets the user of the measurement object.

        Arguments:
                user(str) - - Name of the user.
        """
        self.tags.user = user

    def setProject(self, project):
        """
        Sets the project.

        Arguments:
                project(str) - - Name of the project.
        """
        self.tags.project = project

    def setLogChannels(self, log_channels, override=False):
        """
        Sets the log channels of the measurement.

        Arguments:
                log_channels(list of str or str) - - List of log channel names or a single log channel name. Does not check if the log channels provided actually exist or is valid. Duplicates are removed.
                override(bool) - - If true, overrides the existing log channels. Default is False.
        """

        if isinstance(log_channels, list):

            available_channels = self.getChannelNames()
            for log_channel in log_channels:
                if log_channel not in available_channels:
                    logging.info(
                        log_channel + ' not found in list of channels. Trying to add it.')
                    split_name = log_channel.split(' - ', 1)
                    self.addChannel(split_name[0], split_name[1])

            if override:
                self.log_channels = log_channels
            else:
                current_log_channels = self.getLogChannels()
                # self.scenario['log_channels'] = list(set(current_log_channels) | set(log_channels)) # remove duplicates
                current_log_channels.extend(log_channels)
                self.log_channels = sorted(
                    set(current_log_channels), key=current_log_channels.index)  # remove duplicates
                # self.scenario['log_channels'].extend(log_channels)
        else:
            available_channels = self.getChannelNames()
            if log_channels not in available_channels:
                logging.info(
                    log_channels + ' not found in list of channels. Trying to add it.')
                split_name = log_channels.split(' - ', 1)
                self.addChannel(split_name[0], split_name[1])

            if override:
                self.scenario.log_channels = [log_channels]
            else:
                if log_channels not in self.log_channels:
                    self.log_channels.append(log_channels)

    def updateSettings(self, settings_dict):
        """
        Updates the settings of the measurement scenario.

        Arguments:
        settings_dict(dict) - - dictionary, which keys are the settings names and values are the settings values.

        Available settings and their defaults are:

        'Update instruments at start even if values are unchanged': True,
        'Send values in parallel': True,
        'Limit hardware looping to first step item': False,
        'Number of step items in hardware loop': 1,
        'Only send signal if source instrument has been updated': True,
        'Data compression': 4,
        'Method': 'Nelder-Mead',
        'Max evaluations': 200,
        'Minimization function': 'min(y[0])',
        'Target value': -inf,
        'Relative tolerance': inf,
        'opt-Bayesian-Gaussian-Process: Acquisition function': 'gp_hedge',
        'opt-Bayesian-Gaussian-Process: kappa': 1.96,
        'opt-Bayesian-Gaussian-Process: xi': 0.1,
        """

        for key, value in settings_dict.items():
            if hasattr(self.settings, key):
                setattr(self.settings, key, value)
            else:
                logging.warning(f'Setting {key} not available.')

    def activate_hardware_loop(self, trigger_channel, hardware_loop_depth=1):
        """ Activates hardware looping.

        Args:
            trigger_channel (str): Name of the channel that is used for triggering the hardware loop.
            hardware_loop_depth (int): Limits the hardware looping to n first steps. 0 means no limitiations. Negative values indicate that hardware looping should be disabled.

        """
        if hardware_loop_depth < 0:
            self.updateSettings({
                'hardware_loop': False
            })
            return
        else:
            self.updateSettings({
                'hardware_loop': True
            })

        if hardware_loop_depth < 1:
            self.updateSettings({
                'limit_hardware_looping': False
            })
        else:
            self.updateSettings({
                'limit_hardware_looping': True
            })

        self.updateSettings({
            'n_items_hardware_loop': hardware_loop_depth,
            'arm_trig_mode': True,
            'trig_channel': trigger_channel,
        })

        # XXX
        pass

    def updateInstrumentValueFullName(self, full_name, value):
        """
        This is a wrapper to updateInstrumentValue, which allows you to add full_name.

        Args:
           full_name (str): Full name of the instrument value to be set. Format: "[Instrument name] - [Channel name]"
           value (float): The value of the parameter.
        """
        return self.updateInstrumentValue(None, None, value=value, full_name=full_name)

    def updateInstrumentValuesByDict(self, instrument_values):
        """
            Updates instrument values that are stored in a dictionary.

            The dictionary should have the format:
                instrument_values = {instrument_name: {
                    parameter_name: parameter_value}}
        """
        if not isinstance(instrument_values, dict):
            raise Exception(
                "Instrument value dictionary does not have the correct shape.")
        for instrument_name, parameter_dict in instrument_values.items():
            if not isinstance(parameter_dict, dict):
                raise Exception(
                    "Instrument value dictionary does not have the correct shape.")
            for parameter_name, value in parameter_dict.items():
                self.updateInstrumentValue(
                    instrument_name, parameter_name, value)

    def update_step_values(self,
                           channel_name,
                           values=None,
                           center=None,
                           span=None,
                           start=None,
                           stop=None,
                           n_pts=None,
                           step=None,
                           channel_position=None,
                           range_index=None,
                           override_existing=True
    ):
        """ Updates a step channel, or creates a new channel if it does not exist.

        This function is a more user-friendly version of updateValue. If step channel does not exist, a new channel is created. Calling this function disables equation from the channel. There are multiple ways to provide the values of the step channel. If values is not None, it is used as the source. Second option is center - span. Lowest priority is given to start-stop. The range is created either by using n_pts or a fixed step. If channel_position is given, the step channel is moved to that position in the list of step_channels. If range_index is given, the new values are added to the list of range_items for that step in the given position. If override_existing is True, range_item in that position is overridden, or if range_item is None, all the range_items are overridden.

        Args:
            channel_name (str):
            values (float or list):
            channel_position (int or None):
            range_index (int or None):
            override_existing (bool): If True, all existing range items will be overridden.
        Returns:
            List of step values or None if channel does not exists and cannot be created.
        """
        if values is not None:
            if isinstance(values, Iterable):
                for ii, value in enumerate(values):
                    self.updateValue(channel_name, values,
                                     itemType='SINGLE', step_index=ii)
        logging.warning('update_step_values is not implemented')
        return self.getStepChannelValues(channel_name)

    def updateStepChannelsByDict(self, step_channels, skip_nones=False):
        """
        Updates step channels that are stored in a dictionary.

        The dictionary should have the format:
            step_channels = {channel_name: channel_value}
        Channel value can either be a single value or a dictionary of itemType - value pairs or a list of single values/dictionaries of itemType - value pairs.

        Args:
            step_channels (dict): dictionary if step channels.
            skip_nones (bool): If True, None entries are skipped.
        """
#        for name, value in step_channels.items():
#            #try:
#            if isinstance(value,dict):
#                for dict_key, dict_value in value.items():
#                    self.updateValue(name,dict_value,itemType=dict_key)
#            else:
#                self.updateValue(name, value, itemType='SINGLE')

        for name, value in step_channels.items():
            if isinstance(value, dict):  # Only a single dictionary
                for dict_key, dict_value in value.items():
                    if not (skip_nones and dict_value is None):
                        self.updateValue(name, dict_value, itemType=dict_key)
            elif not isinstance(value, Iterable):
                if not (skip_nones and value is None):
                    self.updateValue(name, value, itemType='SINGLE')
            else:  # list of step items
                # loop over all the steps
                for ii, list_value in enumerate(value):
                    # this step is a dictionary describing the step
                    if isinstance(list_value, dict):
                        for dict_key, dict_value in list_value.items():
                            if not (skip_nones and dict_value is None):
                                self.updateValue(name, dict_value,
                                                 itemType=dict_key, step_index=ii)
                    else:  # this step is just a single point
                        if not (skip_nones and list_value is None):
                            self.updateValue(name, list_value,
                                             itemType='SINGLE', step_index=ii)

    def getOutputPathOfPreviousMeasurement(self):
        """
        Returns the full path of the log file created by the previous measurement.

        Returns None if no previous measurement has been performed.
        """
        print('Warning! getOutputPathOfPreviousMeasurement() is not implmented correctly at the moment.')
        return self.file_out + '.hdf5'

    def update_optimizer_settings(self, optimizer_parameters):
        """ Updates general optimizer parmaters.

        Args:
            optimizer_parameters (dict): dict of setting/value pairs.
        """
        optimizer = self.optimizer  #Labber.config.scenario.Optimizer(**optimizer_parameters)
        for key, value in optimizer_parameters.items():
            setattr(optimizer, key, value)

    def useLabberOptimizer(self, opt_channels, general_settings={}):
        """Activate Labber optimization.
        
        Args:
            opt_channels (dict): Keys are optimization channels, values are optimizer config:
            Example: opt_channels=
                {'param1': {'enabled': True,
                            'initial_step_size': 0.05,
                            'max_value': 1.0,
                            'min_value': 0.0,
                            'precision': 0.001,
                            'start_value': 0.5},
                {'param2': {'enabled': True,
                            'initial_step_size': 10e6,
                            'max_value': 5e9,
                            'min_value': 4e9,
                            'precision': 1e6,
                            'start_value': 4.5e9}}
            general_settings (dict): General optimizer settings, should contain 
            {'max_evaluations': 200,
                'method': 'Nelder-Mead',
                'minimization_function': '1-y[0]',
                'target_value': 0.0,
                'relative_tolerance': 0.02}
        """

        step_channels = self.getStepChannels()

        #go through optimization channels
        for channel, opt_settings in opt_channels.items():
            #cleanup old stuff
            self.removeStepChannel(channel)
            self.updateStepChannelsByDict({channel: 0})

            #go through list of step channels and write settings to Labber optimizer config
            for i, step_ch in enumerate(step_channels):
                if step_ch.channel_name == channel:
                    optimizer = Labber.config.step.OptimizerItem(**opt_settings)
                    step_channels[i].optimizer_config = optimizer

        #general optimizer settings
        self.update_optimizer_settings(general_settings)


def get_full_output_path(output_directory_root, output_file_name):
    """
    Returns the full output path of the next measurement.

    The output path returned is the one that is provided to labber API. Labber API might modify the output path. Therefore this function cannot be used to reliably get the path of previous measurements.
    """
    now = datetime.datetime.now()
    output_directory = output_directory_root + \
        '/{:d}/{:02d}/Data_{:02d}{:02d}/'.format(
            now.year, now.month, now.month, now.day)
    output_path = os.path.join(output_directory, output_file_name)
    return output_path + '.hdf5'
