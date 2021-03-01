from collections import Iterable
from shutil import copyfile
import datetime
import logging
import numpy as np
import sys
import os
import copy
sys.path.append(os.environ.get(
    'LABBERPATH', r'C:\Program Files\Labber\Script'))  # NOQA: E402
import Labber as Labber
from Labber import ScriptTools


class ScriptObject(ScriptTools.MeasurementObject):
    """
    A Labber measurement object based on the JSON dictionary. Contains some additional features.
    """

    def __init__(self, sCfgFileIn, sCfgFileOut, createLog=False):

        self.file_in = sCfgFileIn
        # if sCfgFileOut = None:
        #    output_directory =  +'..\\{:d}\\{:02d}\\Data_{:02d}{:02d}\\'.format(now.year,now.month,now.month,now.day)
        #    output_path = os.path.join(output_directory,output_file_name)
        # else:
        self.file_out = sCfgFileOut
        self.scenario = ScriptTools.load_scenario_as_dict(self.file_in)
        self.master_channel = None

    def performMeasurement(self, return_data=True):
        temp_file_name = os.path.splitext(self.file_out)[0]+'_tmp.labber'
        file_dir_path = os.path.dirname(self.file_out)

        # Create subdirectories if they don't exist.
        if file_dir_path.strip():
            os.makedirs(file_dir_path, exist_ok=True)

        ScriptTools.save_scenario_as_binary(self.scenario, temp_file_name)
        labber_meas_object = ScriptTools.MeasurementObject(
            temp_file_name, self.file_out)
        if self.master_channel is not None:
            labber_meas_object.setMasterChannel(self.master_channel)
        return labber_meas_object.performMeasurement(return_data)

    def save_as_binary(self, filename):
        """ Saves the scipt object as labber binary (.labber file).

        Args:
            filename (str): name of the saved binary.
        """
        ScriptTools.save_scenario_as_binary(self.scenario, filename)

    def setMasterChannel(self, channel_name):
        """
            Sets the master channel.

            When other parameters are updated, creates a lookup table for them with respect to the master channel.

            Does not work correctly at the moment.
        """

        logging.warning('setMasterChannel() is not supported')
        # raise Exception('setMasterChannel not implementd.')
        self.master_channel = channel_name
        # temp_file_name = os.path.splitext(self.file_out)[0]+'_tmp.labber'
        # ScriptTools.save_scenario_as_binary(self.scenario, temp_file_name)
        # labber_meas_object = ScriptTools.MeasurementObject(temp_file_name,self.file_out)
        # labber_meas_object.setMasterChannel(channel_name)

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

        for current_instrument in self.scenario['instruments']:
            if instrument_name == current_instrument['com_config']['name']:
                current_instrument['values'][parameter_name] = value
                return

        logging.warning('Instrument ' + instrument_name + ' not found.')
        return

    def copy_instrument_values(self, source_instrument_name, target_instrument_name):
        """ Copies the instrument values, except for communication related parameters. Instruments have to be of same type.
        """
        source = self.get_instrument(source_instrument_name)
        target = self.get_instrument(target_instrument_name)
        if source is not None and target is not None:
            target['values'] = copy.deepcopy(source['values'])

    def get_instrument(self, instrument_name):
        """ Returns the instrument with given name.

        """
        for current_instrument in self.scenario['instruments']:
            if instrument_name == current_instrument['com_config']['name']:
                return current_instrument
        return None


    def removeStepChannel(self, channel):
        """[summary]
        
        Args:
            channel ([type]): [description]
        """
        for i, c in enumerate(self.scenario['step_channels']):
            if c['channel_name']==channel:
                del self.scenario['step_channels'][i]

    def getStepChannels(self):
        """
        Returns all the step channels (right side of the measurement program).
        """
        return self.scenario['step_channels']

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
            channel_name = channel['channel_name']

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
        return self.scenario['channels']

    def getChannel(self, name):
        """
        Get a reference to channel by its name.

        Arguments:
        name (str) -- Name of the channel.

        Returns:
        Reference to the channel or None if not found.
        """

        # if name in self.getChannelNames():
        for channel in self.getChannels():
            if channel['instrument'] + ' - ' + channel['quantity'] == name:
                return channel
            elif 'name' in channel.keys() and channel['name'] == name:
                return channel
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
        channel_names = []
        for channel in self.getChannels():
            if 'name' in channel.keys():
                channel_names.append(channel['name'])
            else:
                channel_names.append(
                    channel['instrument'] + ' - ' + channel['quantity'])
        return channel_names

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

        instruments = self.scenario['instruments']
        # for current_instrument in self.scenario['instruments']:
        instrument_value_names = []
        for instrument in instruments:
            current_instrument_name = instrument['com_config']['name']
            for value_name in instrument['values'].keys():
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
        log_channels = self.scenario['log_channels']
        return log_channels

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
                instrument_name = channel['instrument']
                quantity_name = channel['quantity']
            else:
                return None
        for current_instrument in self.scenario['instruments']:
            if instrument_name == current_instrument['com_config']['name']:
                return current_instrument['values'][quantity_name]

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
            channel_name = channel['channel_name']

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

        step_channels = self.scenario['step_channels']

        if index < 0:
            index = len(step_channels) + index + 1  # Note the additional +1. -1 moves to the end.
        for ii in range(len(step_channels)):
            channel = step_channels[ii]
            if channel_name == channel['channel_name']:
                step_channels.insert(index, channel)
                if(ii < index):
                    step_channels.pop(ii)
                else:
                    step_channels.pop(ii+1)
                return

    def move_log_channel_to(self, channel_name, index):
        """
        Moves a log channel to the given index in the list of log channels.

        Arguments:
            channel_name(str) - - Name of the log channel
            index(int) - - new position of the channel in the list of channels.
        """

        log_channels = self.scenario['log_channels']

        for i, log_channel in enumerate(log_channels):
            if channel_name==log_channel:
                log_channels.insert(index, log_channel)
                if i<index:
                    log_channels.pop(i)
                else:
                    log_channels.pop(i+1)

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


    def write_lookup(self, lookup_table, channel_name, ref_channel_name, add_current_value=True):
        """ A wrapper to write a lookup table for a channel in equation. Not that this will override the existing equation. The full functionality can be achieved with add_equation() function.

    Args:
        lookup_table (list): Lookup table as list in Labber convention.
        channel_name (Str): Channel name of channel that obtains lookup table.
        ref_channel_name (Str): Reference channel name.
        add_current_value (bool, optional): If True, current value of the channel is added to the lookup table.
    """
        
        if add_current_value:
            self.updateStepChannelsByDict({
                channel_name: {
                    'EQ': 'x+p1',
                    'VARS': {'x': 'Step values', 'p1': {'channel_name': ref_channel_name, 'lookup_x': lookup_table[0], 'lookup_y': lookup_table[1]}}
                }
            })
        else:
            self.updateStepChannelsByDict({channel_name: {'EQ': 'p1', 'VARS': {'x': 'Step values',
                'p1': {'channel_name': ref_channel_name, 'lookup_x': lookup_table[0], 'lookup_y': lookup_table[1]}}}})

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
            if channel is not None and channel['sweep_mode'].lower() != 'Off':
                for step_item in channel['step_items']:
                    step_item['sweep_rate'] = sweep_rate


    def updateValue(self, channel_name, value, itemType='SINGLE', step_index=0):
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

            step_index(int, optional) - - index in the step items which value will be updated.
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
                    channel[param_name] = param_value
                return
            if(itemType == 'EQ'):
                channel['equation'] = value
                set_sweep_parameter(channel, {
                                    'use_relations': True}, 'PARAM', index)
                return
            if(itemType == 'VARS'):
                channel['relation_parameters'] = []
                for dict_param, dict_value in value.items():
                    if(isinstance(dict_value, dict)):
                        if 'interp' in dict_value.keys():
                            interp = dict_value['interp']
                        else:
                            interp = 'Linear'
                            channel['relation_parameters'].append({
                                'variable': dict_param,
                                'channel_name': dict_value['channel_name'],
                                'use_lookup': True,
                                'lookup': {
                                    'interp': interp,
                                    'xdata': np.array(dict_value['lookup_x']),
                                    'ydata': np.array(dict_value['lookup_y'])
                                }
                            })
                    else:
                        channel['relation_parameters'].append(
                            {'variable': dict_param, 'channel_name': dict_value, 'use_lookup': False, 'lookup': None})
                set_sweep_parameter(channel, {
                    'use_relations': True}, 'PARAM', index)

                return

            if(index >= len(channel['step_items'])):
                for ii in range(index-len(channel['step_items'])+1):
                    channel['step_items'].append({
                        'range_type': 'Start - Stop',
                        'step_type': 'Fixed # of pts',
                        'single': 0,
                        'start': 0,
                        'stop': 1,
                        'center': 0.5,
                        'span': 1,
                        'step': 0.1,
                        'n_pts': 1,
                        'interp': 'Linear',
                        'sweep_rate': 0.0
                    })
            if(itemType == 'SINGLE'):
                channel['step_items'][index]['range_type'] = 'Single'
                channel['step_items'][index]['single'] = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'START'):
                channel['step_items'][index]['range_type'] = 'Start - Stop'
                channel['step_items'][index]['start'] = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'STOP'):
                channel['step_items'][index]['range_type'] = 'Start - Stop'
                channel['step_items'][index]['stop'] = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)

            if(itemType == 'SPAN'):
                channel['step_items'][index]['range_type'] = 'Center - Span'
                channel['step_items'][index]['span'] = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'CENTER'):
                channel['step_items'][index]['range_type'] = 'Center - Span'
                channel['step_items'][index]['center'] = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'N_PTS'):
                channel['step_items'][index]['step_type'] = 'Fixed # of pts'
                channel['step_items'][index]['n_pts'] = int(value)
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)
            if(itemType == 'STEP'):
                channel['step_items'][index]['step_type'] = 'Fixed step'
                channel['step_items'][index]['step'] = value
                set_sweep_parameter(channel, {
                    'use_relations': False}, 'PARAM', index)

        step_channels = self.getStepChannels()
        for current_channel in step_channels:
            # print(current_channel['channel_name'],channel_name)
            if(channel_name == current_channel['channel_name']):
                # und match: ' + channel_name + ' and ' + current_channel['channel_name'])
                set_sweep_parameter(current_channel, value,
                                    itemType, step_index)
                # Extend here with advanced options!
                return
        # channel not found, add a new channel.
        new_channel = {'channel_name': channel_name,
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
                       'alternate_direction': False,
                       'step_items': [{
                           'range_type': 'Single',
                           'step_type': 'Fixed # of pts',
                           'single': 0,
                           'start': 0,
                           'stop': 1,
                           'center': 0.5,
                           'span': 1,
                           'step': 0.1,
                           'n_pts': 11,
                           'interp': 'Linear',
                           'sweep_rate': 0.0
                       }],
                       'relation_parameters': [{
                           'variable': 'x',
                           'channel_name': 'Step values',
                           'use_lookup': False, 'lookup': None}],
                       'optimizer_config': {
                           'Enabled': False,
                           'Start value': 7660000000.0,
                           'Initial step size': 1532000000.0,
                           'Min value': 7660000000.0,
                           'Max value': 15320000000.0,
                           'Precision': 766000.0}
                       }
        step_channels.append(new_channel)
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

        target_channel['signal_source'] = source

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
        new_channel = {'channel_name': step_channel_name,
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
                       'alternate_direction': False,
                       'step_items': [{
                           'range_type': 'Single',
                           'step_type': 'Fixed # of pts',
                           'single': value,
                           'start': 0,
                           'stop': 1,
                           'center': 0.5,
                           'span': 1,
                           'step': 0.1,
                           'n_pts': 11,
                           'interp': 'Linear',
                           'sweep_rate': 0.0
                       }],
                       'relation_parameters': [{
                           'variable': 'x',
                           'channel_name': 'Step values',
                           'use_lookup': False, 'lookup': None}],
                       'optimizer_config': {
                           'Enabled': False,
                           'Start value': 7660000000.0,
                           'Initial step size': 1532000000.0,
                           'Min value': 7660000000.0,
                           'Max value': 15320000000.0,
                           'Precision': 766000.0}
                       }

        if step_channel_name in self.getStepChannelNames():
            return
        step_channels = self.getStepChannels()
        step_channels.append(new_channel)

    def copy_step_channel(self, old_channel_name, new_channel_name):
        old_step_channel = self.getStepChannel(old_channel_name)
        if old_step_channel is not None:
            new_step_channel = copy.deepcopy(old_step_channel)
            new_step_channel['channel_name'] = new_channel_name
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

        for channel in self.scenario['channels']:
            if(channel['instrument'] == instrument_name and channel['quantity'] == quantity):
                # Channel already exitst
                if name is not None:
                    self.nameChannel(instrument_name + ' - ' + quantity, name)
                if physical_unit is not None:
                    channel['unit_physical'] = physical_unit
                if instrument_unit is not None:
                    channel['unit_instrument'] = instrument_unit
                return

        channel = {
            'instrument': instrument_name,
            'quantity': quantity,
        }
        if physical_unit is not None:
            channel['unit_physical'] = physical_unit
        if instrument_unit is not None:
            channel['unit_instrument'] = instrument_unit
        if name is not None:
            channel['name'] = name
        self.scenario['channels'].append(channel)

    def nameChannel(self, old_name, new_name):
        """ Gives a nick name to the channel.

        In addition to naming the channel, renames all the step channels and log channels to which
        old_name of the channel refers to.
        Args:
            old_name (str): Full name of the channel in format 'instrument - quantity' or the nick name of the channel.
            nick_name (str): The name with which the channel should be called. If None, nick Name is removed.
        """
        channel = self.getChannel(old_name)

        if channel is not None:
            full_channel_name = channel['instrument'] + ' - ' + channel['quantity']
            if 'name' in channel:
                nick_name = channel['name']
            else:
                nick_name = None
            if new_name is not None:
                channel['name'] = new_name
            else:
                if 'name' in channel:
                    del(channel['name'])
            step_channel = self.getStepChannel(old_name)
            if step_channel is None:
                step_channel = self.getStepChannel(nick_name)
            if step_channel is not None:
                if new_name is not None:
                    step_channel['channel_name'] = new_name
                else:
                    step_channel['channel_name'] = full_channel_name
            for ii, log_channel in enumerate(self.getLogChannels()):
                if log_channel == old_name or log_channel == nick_name:
                    if new_name is not None:
                        self.getLogChannels()[ii] = new_name
                    else:
                        self.getLogChannels()[ii] = channel['instrument'] + ' - ' + channel['quantity']

    def setParameter(self, name, value):
        """
        Sets a parameter value pair for the measurement object.

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
        if name in self.scenario['parameters']:
            self.scenario['parameters'][name] = value
        else:
            logging.warning('Parameter ' + name + ' does not exist.')

    def setComment(self, comment):
        """
        Add a comment to the measurement.
        """
        if comment is not None:
            self.scenario['parameters']['comment'] = comment

    def setTags(self, tags, override=True):
        """
        Sets the tags for the measurement object.

        Arguments:
                tags(list of str or str) - - List tags to be set or a single tag
                override(bool) - - If True, all the existing tags will be overridden. If False, appends tags.
        """
        if isinstance(tags, list):

            if override:
                self.scenario['tags']['tags'] = tags
            else:
                self.scenario['tags']['tags'].extend(tags)
        else:
            if override:
                self.scenario['tags']['tags'] = [tags]
            else:
                self.scenario['tags']['tags'].append(tags)

    def setUser(self, user):
        """
        Sets the user of the measurement object.

        Arguments:
                user(str) - - Name of the user.
        """
        self.scenario['tags']['user'] = user

    def setProject(self, project):
        """
        Sets the project.

        Arguments:
                project(str) - - Name of the project.
        """
        self.scenario['tags']['project'] = project

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
                self.scenario['log_channels'] = log_channels
            else:
                current_log_channels = self.getLogChannels()
                # self.scenario['log_channels'] = list(set(current_log_channels) | set(log_channels)) # remove duplicates
                current_log_channels.extend(log_channels)
                self.scenario['log_channels'] = sorted(
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
                self.scenario['log_channels'] = [log_channels]
            else:
                if log_channels not in self.scenario['log_channels']:
                    self.scenario['log_channels'].append(log_channels)

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

        self.scenario['settings'].update(settings_dict)

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

    def updateValuesByArrays(self, globals_path, tags=None, comment=None, user=None, project=None, globals_file_name=globals, local_step_channels={}, looped_variables={}, local_parameters={}, local_instrument_values={}, log_channels=None, override_log_channels=False,local_settings={}):
        """
        Updates the measurement object using parameter dictionaries.
        """
        sys.path.append(globals_path)
        gl = __import__(globals_file_name)

        self.updateInstrumentValuesByDict(gl.instrument_values)
        self.updateInstrumentValuesByDict(local_instrument_values)

        self.updateStepChannelsByDict(gl.step_channels)
        self.updateStepChannelsByDict(local_step_channels)
        self.updateSettings(local_settings)
        self.updateStepChannelsByDict(looped_variables)
#        for instrument_name,parameter_dict in gl.instrument_values.items():
#            for parameter_name,value in parameter_dict.items():
#                self.updateInstrumentValue(instrument_name,parameter_name,value)
#
#        for instrument_name,parameter_dict in local_instrument_values.items():
#            for parameter_name,value in parameter_dict.items():
#                self.updateInstrumentValue(instrument_name,parameter_name,value)


#        for name, value in gl.step_channels.items():
#            #try:
#            if isinstance(value,dict):
#                for dict_key, dict_value in value.items():
#                    self.updateValue(name,dict_value,itemType=dict_key)
#            else:
#                self.updateValue(name, value, itemType='SINGLE')
#
#
#            #measurement_object.updateValue(name, value, itemType='SINGLE')
#            #except Exception as e:
#            #   logging.warning('Unable to set variable ' + name +'.') # Does nothing at the moment
#
#
#        for name, value in local_step_channels.items():
#            if isinstance(value,dict):
#                for dict_key, dict_value in value.items():
#                    self.updateValue(name,dict_value,itemType=dict_key)
#            else:
#                self.updateValue(name, value, itemType='SINGLE')

#        for name,value in looped_variables.items():
#            if isinstance(value,dict): # Only a single dictionary
#                for dict_key, dict_value in value.items():
#                    self.updateValue(name,dict_value,itemType=dict_key)
#            else: # list of step items
#                for ii,list_value in enumerate(value): # loop over all the steps
#                    if isinstance(list_value,dict): # this step is a dictionary describing the step
#                        for dict_key, dict_value in list_value.items():
#                            self.updateValue(name,dict_value,itemType=dict_key,step_index=ii)
#                    else: # this step is just a single point
#                        self.updateValue(name,list_value,itemType='SINGLE',step_index=ii)

        counter = 0
        for name in looped_variables.keys():
            self.move_channel_to(name, counter)
            counter = counter+1

        if gl.log_channels is not None:
            self.setLogChannels(log_channels, override_log_channels)

        if log_channels is not None:
            self.setLogChannels(log_channels, override_log_channels)

        for name, value in gl.parameters.items():
            self.setParameter(name, value)

        for name, value in local_parameters.items():
            self.setParameter(name, value)

        if user is not None:
            self.setUser(user)
        if comment is not None:
            self.setComment(comment)
        if project is not None:
            self.setProject(project)
        if tags is not None:
            self.setTags(tags)

    def getOutputPathOfPreviousMeasurement(self):
        """
        Returns the full path of the log file created by the previous measurement.

        Returns None if no previous measurement has been performed.
        """
        print('Warning! getOutputPathOfPreviousMeasurement() is not implmented correctly at the moment.')
        return self.file_out + '.hdf5'


    def useLabberOptimizer(self, opt_channels, general_settings):
        """Activate Labber optimization.
        
        Args:
            opt_channels (dict): Keys are optimization channels, values are optimizer config:
            Example: opt_channels=
                {'param1': {'Enabled': True,
                            'Initial step size': 0.05,
                            'Max value': 1.0,
                            'Min value': 0.0,
                            'Precision': 0.001,
                            'Start value': 0.5},
                {'param2': {'Enabled': True,
                            'Initial step size': 10e6,
                            'Max value': 5e9,
                            'Min value': 4e9,
                            'Precision': 1e6,
                            'Start value': 4.5e9}}
            general_settings (dict): General optimizer settings, should contain 
            {'Max evaluations': 200,
                'Method': 'Nelder-Mead',
                'Minimization function': '1-y[0]',
                'Target value': 0.0,
                'Relative tolerance': 0.02}
        """

        step_channels=self.getStepChannels()

        #go through optimization channels
        for channel, opt_settings in opt_channels.items():
            #cleanup old stuff
            self.removeStepChannel(channel)
            self.updateStepChannelsByDict({channel: 0})

            #go through list of step channels and write settings to Labber optimizer config
            for i, step_ch in enumerate(step_channels):
                if step_ch['channel_name']==channel:
                    step_channels[i]['optimizer_config'].update(opt_settings)

        #general optimizer settings
        self.scenario['settings'].update(general_settings)


def updateAndPerformMeasurement_old(template_path, output_directory_root, output_file_name, tags, comment, globals_path, globals_file_name=globals, local_step_channels={}, looped_variables={}, local_parameters={}, local_instrument_values={}, log_channels=None):
    now = datetime.datetime.now()
    output_directory = os.path.join(output_directory_root, '/{:d}/{:02d}/Data_{:02d}{:02d}/'.format(
        now.year, now.month, now.month, now.day))
    output_path = os.path.join(output_directory, output_file_name)
    output_log_file = os.path.join(output_directory, output_file_name) + '.log'
    print(output_log_file)
    script_name = os.path.realpath(__file__)

    sys.path.append(globals_path)
    gl = __import__(globals_file_name)

    now = datetime.datetime.now()
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    logging.basicConfig(filename=output_log_file,
                        format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)
    root = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    root.addHandler(handler)
    # measurement_object=ScriptTools.load_scenario_as_dict(template_path)
    measurement_object = ScriptObject(template_path, output_path)

    logging.info("Starting measurement " + output_file_name)

    for instrument_name, parameter_dict in gl.instrument_values.items():
        for parameter_name, value in parameter_dict.items():
            measurement_object.updateInstrumentValue(
                instrument_name, parameter_name, value)

    for instrument_name, parameter_dict in local_instrument_values.items():
        for parameter_name, value in parameter_dict.items():
            measurement_object.updateInstrumentValue(
                instrument_name, parameter_name, value)

    for name, value in gl.step_channels.items():
        # try:
        if isinstance(value, dict):
            for dict_key, dict_value in value.items():
                measurement_object.updateValue(
                    name, dict_value, itemType=dict_key)
        else:
            measurement_object.updateValue(name, value, itemType='SINGLE')
        # measurement_object.updateValue(name, value, itemType='SINGLE')
        # except Exception as e:
        #   logging.warning('Unable to set variable ' + name +'.') # Does nothing at the moment

    for name, value in local_step_channels.items():
        if isinstance(value, dict):
            for dict_key, dict_value in value.items():
                measurement_object.updateValue(
                    name, dict_value, itemType=dict_key)
        else:
            measurement_object.updateValue(name, value, itemType='SINGLE')

    for name, value in looped_variables.items():
        if isinstance(value, dict):  # Only a single dictionary
            for dict_key, dict_value in value.items():
                measurement_object.updateValue(
                    name, dict_value, itemType=dict_key)
        else:  # list of step items
            for ii, list_value in enumerate(value):  # loop over all the steps
                # this step is a dictionary describing the step
                if isinstance(list_value, dict):
                    for dict_key, dict_value in list_value.items():
                        measurement_object.updateValue(
                            name, dict_value, itemType=dict_key, step_index=ii)
                else:  # this step is just a single point
                    measurement_object.updateValue(
                        name, list_value, itemType='SINGLE', step_index=ii)

    counter = 0
    for name in looped_variables.keys():
        measurement_object.move_channel_to(name, counter)
        counter = counter+1

    if gl.log_channels is not None:
        measurement_object.setLogChannels(log_channels)

    if log_channels is not None:
        measurement_object.setLogChannels(log_channels)

    for name, value in gl.parameters.items():
        measurement_object.setParameter(name, value)

    for name, value in local_parameters.items():
        measurement_object.setParameter(name, value)

    # sys.exit("Error message")
    # Copy the currently running script
    copyfile(script_name, output_path + ".py")
    # raise Exception('Debugging')
    measurement_object.performMeasurement()

    Lfile = Labber.LogFile(output_path + '.hdf5')

    tags.extend(gl.tags)
    Lfile.setTags(tags)
    Lfile.setUser(gl.user)
    Lfile.setComment(comment)
    Lfile.setProject(gl.project)


def updateAndPerformMeasurement(template_path, output_directory_root, output_file_name, tags, comment, globals_path, globals_file_name=globals, local_step_channels={}, looped_variables={}, local_parameters={}, local_instrument_values={}, log_channels=None,local_settings={}):
    print('Entering')
    now = datetime.datetime.now()
    output_directory = output_directory_root + \
        '/{:d}/{:02d}/Data_{:02d}{:02d}/'.format(
            now.year, now.month, now.month, now.day)
    output_path = os.path.join(output_directory, output_file_name)
    output_log_file = os.path.join(output_directory, output_file_name) + '.log'
    print(output_log_file)
    print('Output path: ' + output_path)
    script_name = os.path.realpath(__file__)

    sys.path.append(globals_path)
    gl = __import__(globals_file_name)

    now = datetime.datetime.now()
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    logging.basicConfig(filename=output_log_file,
                        format='%(asctime)s:%(levelname)s:%(message)s', level=logging.DEBUG)
    root = logging.getLogger()
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    root.addHandler(handler)
    # measurement_object=ScriptTools.load_scenario_as_dict(template_path)

    logging.info("Starting measurement " + output_file_name)

    tags.extend(gl.tags)
    measurement_object = ScriptObject(template_path, output_path)
    measurement_object.updateValuesByArrays(globals_path, tags=tags, comment=comment, user=gl.user, project=gl.project, globals_file_name=globals_file_name, local_step_channels=local_step_channels,
                                            looped_variables=looped_variables, local_parameters=local_parameters, local_instrument_values=local_instrument_values, log_channels=log_channels, override_log_channels=False,local_settings=local_settings)
    return measurement_object.performMeasurement()


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
