"""
Testing Module for config_utils.py

This testing module is dedicated to thoroughly testing the core functions
provided by the config_utils.py module.  It focuses on ensuring the accuracy of
file loading, the correctness of configuration file updates, logical
functionality tests, and more.

The primary goals of these tests are to validate the proper handling of
configuration files, the expected behavior of configuration updates, and the
correctness of the underlying logic in the config_utils module.
"""


import pytest

from cytosnake.utils.config_utils import load_configs, update_config


@pytest.mark.positive
def test_add_new_entry(config_file_test):
    """Postive test for entering a new config entry into a config file"""
    # set new key and value
    new_key = "test_key"
    new_value = "new_value"

    # update configs
    update_config(config_file_test, new_key=new_key, new_value=new_value)

    # load config
    configs = load_configs(config_file_test)

    # test
    assert new_key in list(configs.keys())
    assert configs[new_key] == new_value


@pytest.mark.negative
def test_overwriting_key():
    """Overwrite values from existing keys"""

    # set new key and value
    new_key = "test_key"
    new_value = "new_value"

    # update configs
    update_config(config_file_test, new_key=new_key, new_value=new_value)

    # load config
    configs = load_configs(config_file_test)

    # test
    assert new_key in list(configs.keys())
    assert configs[new_key] == new_value
